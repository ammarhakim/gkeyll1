/**
 * @file	LcLagrangeTensorElement.cpp
 *
 * @brief       Lagrange tensor-product element.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLagrangeTensorElement.h>
#include <LcStructuredGridBase.h>

namespace Lucee
{
  template <> const char *LagrangeTensorElement<1>::id = "LagrangeTensor";
  template <> const char *LagrangeTensorElement<2>::id = "LagrangeTensor";
  template <> const char *LagrangeTensorElement<3>::id = "LagrangeTensor";
  template <> const char *LagrangeTensorElement<4>::id = "LagrangeTensor";
  template <> const char *LagrangeTensorElement<5>::id = "LagrangeTensor";

  template <unsigned NDIM>
  LagrangeTensorElement<NDIM>::LagrangeTensorElement()
    : Lucee::NodalFiniteElementIfc<NDIM>(1), nodeSeq(Lucee::Region<NDIM, int>()),
      local2Global(Lucee::Region<NDIM, int>())
  {
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::NodalFiniteElementIfc<NDIM>::readInput(tbl);

// get polynomial order to use
    if (tbl.hasNumber("polyOrder"))
    {
      unsigned polyOrder = tbl.getNumber("polyOrder");
      for (unsigned d=0; d<NDIM; ++d)
        numNodes[d] = ((unsigned) polyOrder) + 1;
    }
    else if (tbl.hasNumVec("polyOrder"))
    {
      std::vector<double> po = tbl.getNumVec("polyOrder");
      if (po.size() != NDIM)
      {
        Lucee::Except lce("LagrangeTensorElement::readInput: 'polyOrder' should have ");
        lce << NDIM << " entries. Provided " << po.size() << " instead.";
        throw lce;
      }
      for (unsigned d=0; d<NDIM; ++d)
        numNodes[d] = ((unsigned) po[d]) + 1;
    }
    else
    {
      Lucee::Except lce("LagrangeTensorElement::readInput: Must provide 'polyOrder' ");
      lce << "as scalar or table with " << NDIM << " entries.";
      throw lce;
    }

// get node location
    nodeLoc = Lucee::UNIFORM; // by default assume uniformly spaced nodes
    if (tbl.hasString("nodeLocation"))
    {
      std::string nl = tbl.getString("nodeLocation");
      if (nl == "lobatto")
        nodeLoc = Lucee::LOBATTO;
      else if (nl == "gaussian")
        nodeLoc = Lucee::GAUSSIAN;
      else if (nl == "uniform")
        nodeLoc = Lucee::UNIFORM;
      else
      {
        Lucee::Except lce("LagrangeTensorElement::readInput: Node location '");
        lce << nl << "' not recognized.";
        throw lce;
      }
    }
// initialize calculator object
    basisCalc.calc(nodeLoc, numNodes);

// initialize number of local nodes
    unsigned nlocal = basisCalc.getNumNodes();
    this->setNumNodes(nlocal);

// store exclusively owned node indices
    basisCalc.getExclusiveNodeIndices(exclusiveNodes);
// store list of indices for exclusively owned nodes
    basisCalc.getExclusiveNodeNdimIndices(exclusiveNodeIndices.indices);

// construct sequencer for use in mapping functions
    int lower[NDIM], upper[NDIM];
    for (unsigned d=0; d<NDIM; ++d) 
    {
      lower[d] = 0;
      upper[d] = numNodes[d];
    }
    Lucee::Region<NDIM, int> nodeRgn(lower, upper);
    nodeSeq = Lucee::RowMajorSequencer<NDIM>(nodeRgn);
    
// store list of indices of nodes on faces
    for (unsigned d=0; d<NDIM; ++d)
    {
      std::vector<int> idx(NDIM);

// lower face indices
      Lucee::RowMajorSequencer<NDIM> lowerSeq(
        nodeRgn.resetBounds(d, nodeRgn.getLower(d), nodeRgn.getLower(0)+1));
      while (lowerSeq.step())
      {
        lowerSeq.fillWithIndex(&idx[0]);
        lowerIndices[d].indices.push_back(idx);
      }

// upper face indices
      Lucee::RowMajorSequencer<NDIM> upperSeq(
        nodeRgn.resetBounds(d, nodeRgn.getUpper(d)-1, nodeRgn.getUpper(0)));
      while (upperSeq.step())
      {
        upperSeq.fillWithIndex(&idx[0]);
        upperIndices[d].indices.push_back(idx);
      }
    }

// get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->template getGrid<Lucee::StructuredGridBase<NDIM> >();
    Lucee::Region<NDIM, int> gridRgn = grid.getGlobalRegion();

// create indexer for mapping local->global nodes
    if (nodeLoc == GAUSSIAN)
    {
      int lower[NDIM], upper[NDIM];
      for (unsigned d=0; d<NDIM; ++d)
      {
        lower[d] = 0;
        upper[d] = numNodes[d]*gridRgn.getShape(d);
      }
      local2GlobalRgn = Lucee::Region<NDIM, int>(lower, upper);
      local2Global = Lucee::RowMajorIndexer<NDIM>(local2GlobalRgn);

// compute strides
      for (unsigned d=0; d<NDIM; ++d)
        lgStrides[d] = numNodes[d];
    }
    else
    {
      int lower[NDIM], upper[NDIM];
      for (unsigned d=0; d<NDIM; ++d)
      {
        lower[d] = 0;
        upper[d] = (numNodes[d]-1)*gridRgn.getShape(d) + 1;
      }
      local2GlobalRgn = Lucee::Region<NDIM, int>(lower, upper);
      local2Global = Lucee::RowMajorIndexer<NDIM>(local2GlobalRgn);

// compute strides
      for (unsigned d=0; d<NDIM; ++d)
        lgStrides[d] = numNodes[d]-1;
    }

    if (nodeLoc == GAUSSIAN)
    {
      numGlobalNodes = 1;
      for (unsigned d=0; d<NDIM; ++d)
        numGlobalNodes *= numNodes[d]*gridRgn.getShape(d); // no nodes on faces
    }
    else
    {
// compute total number of global nodes
      numGlobalNodes = 1;
      for (unsigned d=0; d<NDIM; ++d)
        numGlobalNodes *= 1 + (numNodes[d]-1)*gridRgn.getShape(d);
    }

// scaling factors for use
    double dx[NDIM], vol2 = 1.0;
    for (unsigned d=0; d<NDIM; ++d)
    {
      dx[d] = grid.getDx(d);
      vol2 *= 0.5*dx[d];
    }

// mass matrix
    mass = Lucee::Matrix<double>(nlocal, nlocal);
    basisCalc.getMassMatrix(mass);
    mass *= vol2;

// grad-stiff matrices
    for (unsigned d=0; d<NDIM; ++d)
    {
      gradStiff[d] = Lucee::Matrix<double>(nlocal, nlocal);
      basisCalc.getGradStiffnessMatrix(d, gradStiff[d]);
      gradStiff[d] *= 2.0/dx[d]*vol2;
    }

// face mass-matrices
    for (unsigned d=0; d<NDIM; ++d)
    {
      unsigned nf = basisCalc.getNumSurfLowerNodes(d);
      lowerFace[d] = Lucee::Matrix<double>(nlocal, nf);
      basisCalc.getLowerFaceMassMatrix(d, lowerFace[d]);

      nf = basisCalc.getNumSurfUpperNodes(d);
      upperFace[d] = Lucee::Matrix<double>(nlocal, nf);
      basisCalc.getUpperFaceMassMatrix(d, upperFace[d]);

// compute factor to bring into physical space: this is promotional to
// area of face normal to direction 'd'.
      double fact = 1.0;
      for (unsigned dd=0; dd<NDIM; ++dd)
        if (dd != d)
          fact *= 0.5*dx[dd];

// scale matrices appropriately
      lowerFace[d] *= fact;
      upperFace[d] *= fact;
    }

// stiffness matrix
    stiff = Lucee::Matrix<double>(nlocal, nlocal);
    basisCalc.getStiffnessMatrix(dx, stiff);
    stiff *= vol2;

    double eta[NDIM];
// compute coordinates of nodes relative to lower-left vertex
    localNodeCoords = Lucee::Matrix<double>(nlocal, NDIM);
    localNodeCoords = 0.0;
    for (unsigned i=0; i<nlocal; ++i)
    {
      basisCalc.fillWithNodeCoordinate(i, eta);
      for (unsigned d=0; d<NDIM; ++d)
        localNodeCoords(i,d) = (eta[d]+1)*0.5*dx[d];
    }

// compute quadrature data for volume and surface quadrature
    volumeQuad.interp = Lucee::Matrix<double>(nlocal, nlocal);
    volumeQuad.ordinates = Lucee::Matrix<double>(nlocal, NDIM);
    volumeQuad.weights.resize(nlocal);

    basisCalc.getGaussQuadData(volumeQuad.interp, volumeQuad.ordinates,
      volumeQuad.weights);

// normalize weights
    for (unsigned i=0; i<nlocal; ++i)
      for (unsigned d=0; d<NDIM; ++d)
        volumeQuad.weights[i] *= 0.5*dx[d];

    for (unsigned d=0; d<NDIM; ++d)
    {
// lower surface quadrature
      unsigned nlf = basisCalc.getNumSurfLowerNodes(d);
      lowerSurfQuad[d].interp = Lucee::Matrix<double>(nlf, nlocal);
      lowerSurfQuad[d].ordinates = Lucee::Matrix<double>(nlf, NDIM);
      lowerSurfQuad[d].weights.resize(nlf);
      
      basisCalc.getSurfLowerGaussQuadData(d, lowerSurfQuad[d].interp, lowerSurfQuad[d].ordinates,
        lowerSurfQuad[d].weights);

// normalize weights
      for (unsigned i=0; i<nlf; ++i)
        for (unsigned dd=0; dd<NDIM; ++dd)
          if (dd != d)
            lowerSurfQuad[d].weights[i] *= 0.5*dx[dd];

// upper surface quadrature
      unsigned nuf = basisCalc.getNumSurfUpperNodes(d);
      upperSurfQuad[d].interp = Lucee::Matrix<double>(nuf, nlocal);
      upperSurfQuad[d].ordinates = Lucee::Matrix<double>(nuf, NDIM);
      upperSurfQuad[d].weights.resize(nuf);
      
      basisCalc.getSurfUpperGaussQuadData(d, upperSurfQuad[d].interp, upperSurfQuad[d].ordinates,
        upperSurfQuad[d].weights);

// normalize weights
      for (unsigned i=0; i<nuf; ++i)
        for (unsigned dd=0; dd<NDIM; ++dd)
          if (dd != d)
            upperSurfQuad[d].weights[i] *= 0.5*dx[dd];
    }

// fetch nodal quadrature weights
    nodalWeights.resize(nlocal);
    basisCalc.getWeights(nodalWeights);
// normalize them
    for (unsigned i=0; i<nlocal; ++i)
      nodalWeights[i] *= vol2;
  }

  template <unsigned NDIM>
  unsigned
  LagrangeTensorElement<NDIM>::getNumSurfLowerNodes(unsigned dir) const
  {
    return basisCalc.getNumSurfLowerNodes(dir);
  }

  template <unsigned NDIM>
  unsigned
  LagrangeTensorElement<NDIM>::getNumSurfUpperNodes(unsigned dir) const
  {
    return basisCalc.getNumSurfUpperNodes(dir);
  }

  template <unsigned NDIM>
  unsigned
  LagrangeTensorElement<NDIM>::getNumGlobalNodes() const
  {
    return numGlobalNodes;
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getExclusiveNodeIndices(std::vector<int>& ndIds)
  {
    ndIds = exclusiveNodes;
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getLocalToGlobal(std::vector<int>& lgMap) const
  {
    int nidx[NDIM], idx[NDIM];
    unsigned count = 0;
    nodeSeq.reset();
    while (nodeSeq.step())
    {
      nodeSeq.fillWithIndex(nidx); // local node index
      for (unsigned d=0; d<NDIM; ++d)
        idx[d] = nidx[d] + lgStrides[d]*this->currIdx[d];
      lgMap[count++] = local2Global.getIndex(idx);
    }
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getSurfLowerLocalToGlobal(unsigned dir, std::vector<int>& lgMap) const
  {
    int idx[NDIM];
    unsigned count = 0;
    for (unsigned n=0; n<lowerIndices[dir].indices.size(); ++n)
    {
      for (unsigned d=0; d<NDIM; ++d)
        idx[d] = lgStrides[d]*this->currIdx[d] + lowerIndices[dir].indices[n][d];
      lgMap[count++] = local2Global.getIndex(idx);
    }
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getSurfUpperLocalToGlobal(unsigned dir, std::vector<int>& lgMap) const
  {
    int idx[NDIM];
    unsigned count = 0;
    for (unsigned n=0; n<upperIndices[dir].indices.size(); ++n)
    {
      for (unsigned d=0; d<NDIM; ++d)
        idx[d] = lgStrides[d]*this->currIdx[d] + upperIndices[dir].indices[n][d];
      lgMap[count++] = local2Global.getIndex(idx);
    }
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getSurfLowerNodeNums(unsigned dir, std::vector<int>& nodeNum) const
  {
    basisCalc.getSurfLowerNodeNums(dir, nodeNum);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getSurfUpperNodeNums(unsigned dir, std::vector<int>& nodeNum) const
  {
    basisCalc.getSurfUpperNodeNums(dir, nodeNum);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getNodalCoordinates(Lucee::Matrix<double>& nodeCoords)
  {
// get grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->template getGrid<Lucee::StructuredGridBase<NDIM> >();
// set index and get lower-left vertex coordinates
    grid.setIndex(this->currIdx);
    double xn[7]; // THIS HARD-CODING SEEMS BAD, BUT FOR NOW IT IS (10/26/2012)
    grid.getVertex(xn);
// compute nodal coordinates
    for (unsigned i=0; i<this->getNumNodes(); ++i)
    {
      for (unsigned d=0; d<NDIM; ++d)
        nodeCoords(i,d) = xn[d] + localNodeCoords(i,d);
    }
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getWeights(std::vector<double>& w)
  {
    w.clear(); w.resize(this->template getNumNodes());
    w = nodalWeights;
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getSurfUpperWeights(unsigned dir, std::vector<double>& w)
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getSurfUpperWeights(dir, w);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getSurfLowerWeights(unsigned dir, std::vector<double>& w)
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getSurfLowerWeights(dir, w);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getMassMatrix(Lucee::Matrix<double>& NjNk) const
  {
    NjNk.copy(mass);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getLowerFaceMassMatrix(unsigned dir, Lucee::Matrix<double>& NjNk) const
  {
    NjNk.copy(lowerFace[dir]);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getUpperFaceMassMatrix(unsigned dir, Lucee::Matrix<double>& NjNk) const
  {
    NjNk.copy(upperFace[dir]);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getStiffnessMatrix(Lucee::Matrix<double>& DNjDNk) const
  {
    DNjDNk.copy(stiff);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getGradStiffnessMatrix(unsigned dir, Lucee::Matrix<double>& DNjNk) const
  {
    DNjNk.copy(gradStiff[dir]);
  }

  template <unsigned NDIM>
  unsigned
  LagrangeTensorElement<NDIM>::getNumGaussNodes() const
  {
    return basisCalc.getNumNodes();
  }

  template <unsigned NDIM>
  unsigned
  LagrangeTensorElement<NDIM>::getNumSurfGaussNodes() const
  {
    return basisCalc.getNumSurfUpperNodes(0); // ASSUMPTION: ALL FACES HAVE SAME NUMBER OF SURFACE NODES
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getGaussQuadData(Lucee::Matrix<double>& interpMat,
    Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const
  {
    unsigned nm = NDIM>3 ? 3 : NDIM; // needed for now due to assumption on 3D coordinates
    interpMat.copy(volumeQuad.interp);
    for (unsigned i=0; i<this->getNumGaussNodes(); ++i)
    {
      weights[i] = volumeQuad.weights[i];
      for (unsigned d=0; d<nm; ++d)
        ordinates(i,d) = volumeQuad.ordinates(i,d);
    }
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getSurfLowerGaussQuadData(unsigned dir, Lucee::Matrix<double>& interpMat,
    Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const
  {
    unsigned nm = NDIM>3 ? 3 : NDIM; // needed for now due to assumption on 3D coordinates
    interpMat.copy(lowerSurfQuad[dir].interp);
    for (unsigned i=0; i<this->getNumSurfGaussNodes(); ++i)
    {
      weights[i] = lowerSurfQuad[dir].weights[i];
      for (unsigned d=0; d<nm; ++d)
        ordinates(i,d) = lowerSurfQuad[dir].ordinates(i,d);
    }
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getSurfUpperGaussQuadData(unsigned dir, Lucee::Matrix<double>& interpMat,
    Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const
  {
    unsigned nm = NDIM>3 ? 3 : NDIM; // needed for now due to assumption on 3D coordinates
    interpMat.copy(upperSurfQuad[dir].interp);
    for (unsigned i=0; i<this->getNumSurfGaussNodes(); ++i)
    {
      weights[i] = upperSurfQuad[dir].weights[i];
      for (unsigned d=0; d<nm; ++d)
        ordinates(i,d) = upperSurfQuad[dir].ordinates(i,d);
    }
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getMomentMatrix(unsigned p, Lucee::Matrix<double>& momMat) const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getMomentMatrix(p, momMat);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::extractFromField(const Lucee::Field<NDIM, double>& fld,
    std::vector<double>& data)
  {
// create a region with each side being of size 2
    int lower[NDIM], upper[NDIM];
    for (unsigned d=0; d<NDIM; ++d)
    {
      lower[d] = this->currIdx[d];
      upper[d] = this->currIdx[d]+2;
    }
    Lucee::RowMajorSequencer<NDIM> extRgnSeq(
      Lucee::Region<NDIM, int>(lower, upper));
// create box for mapping to linear index in cell
    for (unsigned d=0; d<NDIM; ++d)
    {
      lower[d] = lgStrides[d]*this->currIdx[d];
      upper[d] = lower[d] + lgStrides[d] + 1;
    }
    Lucee::Region<NDIM, int> extRgn(lower, upper);
    Lucee::RowMajorIndexer<NDIM> extRgnIdx(extRgn);

    Lucee::ConstFieldPtr<double> fldPtr = fld.createConstPtr();
    unsigned nn = exclusiveNodeIndices.indices.size();
    unsigned maxNodes = this->template getNumNodes();
    
    int idx[NDIM], idxN[NDIM];
// loop over region to copy
    while (extRgnSeq.step())
    {
      extRgnSeq.fillWithIndex(idx);
// set pointer to this cell
      fld.setPtr(fldPtr, idx);
// for each exclusively owned node in current cell, check if it should
// be copied and if so, copy it
      for (unsigned n=0; n<nn; ++n)
      {
        for (unsigned d=0; d<NDIM; ++d)
          idxN[d] = lgStrides[d]*idx[d] + exclusiveNodeIndices.indices[n][d];
        if (extRgn.isInside(idxN))
        {
          unsigned loc = extRgnIdx.getIndex(idxN);
          data[loc] = fldPtr[n];
        }
      }
    }
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::copyAllDataFromField(const Lucee::Field<NDIM, double>& fld, double *data)
  {
// region to copy
    Lucee::Region<NDIM, int> rgn =
      this->template getGrid<Lucee::StructuredGridBase<NDIM> >().getLocalRegion();
// extend region by one extra layer of cells on top edges
    int lowerExt[NDIM], upperExt[NDIM];
    for (unsigned d=0; d<NDIM; ++d)
    {
      lowerExt[d] = 0;
      upperExt[d] = 1;
    }
    Lucee::Region<NDIM, int> extRgn = rgn.extend(lowerExt, upperExt);
    Lucee::RowMajorSequencer<NDIM> extRgnSeq(extRgn);

    Lucee::ConstFieldPtr<double> fldPtr = fld.createConstPtr();
    unsigned nn = exclusiveNodeIndices.indices.size();
    unsigned maxNodes = getNumGlobalNodes();
    
    int idx[NDIM], idxN[NDIM];
// loop over region to copy
    while (extRgnSeq.step())
    {
      extRgnSeq.fillWithIndex(idx);
// set pointer to this cell
      fld.setPtr(fldPtr, idx);
// for each exclusively owned node in current cell, check if it should
// be copied and if so, copy it
      for (unsigned n=0; n<nn; ++n)
      {
        for (unsigned d=0; d<NDIM; ++d)
          idxN[d] = lgStrides[d]*idx[d] + exclusiveNodeIndices.indices[n][d];
// check if node belongs to global nodal box
        if (local2GlobalRgn.isInside(idxN))
        {
          unsigned loc = local2Global.getIndex(idxN);
          data[loc] = fldPtr[n];
        }
      }
    }
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::copyAllDataToField(const double *data, Lucee::Field<NDIM, double>& fld)
  {
// region to copy
    Lucee::Region<NDIM, int> rgn =
      this->template getGrid<Lucee::StructuredGridBase<NDIM> >().getLocalRegion();
// extend region by one extra layer of cells on top edges
    int lowerExt[NDIM], upperExt[NDIM];
    for (unsigned d=0; d<NDIM; ++d)
    {
      lowerExt[d] = 0;
      upperExt[d] = 1;
    }
    Lucee::Region<NDIM, int> extRgn = rgn.extend(lowerExt, upperExt);
    Lucee::RowMajorSequencer<NDIM> extRgnSeq(extRgn);

    Lucee::FieldPtr<double> fldPtr = fld.createPtr();
    unsigned nn = exclusiveNodeIndices.indices.size();
    unsigned maxNodes = getNumGlobalNodes();
    
    int idx[NDIM], idxN[NDIM];
// loop over region to copy
    while (extRgnSeq.step())
    {
      extRgnSeq.fillWithIndex(idx);
// set pointer to this cell
      fld.setPtr(fldPtr, idx);
// for each exclusively owned node in current cell, check if it should
// be copied and if so, copy it
      for (unsigned n=0; n<nn; ++n)
      {
        for (unsigned d=0; d<NDIM; ++d)
          idxN[d] = lgStrides[d]*idx[d] + exclusiveNodeIndices.indices[n][d];
// check if node belongs to global nodal box
        if (local2GlobalRgn.isInside(idxN))
        {
          unsigned loc = local2Global.getIndex(idxN);
          fldPtr[n] = data[loc];
        }
      }
    }
  }

// instantiations
  template class LagrangeTensorElement<1>;
  template class LagrangeTensorElement<2>;
  template class LagrangeTensorElement<3>;
  template class LagrangeTensorElement<4>;
  template class LagrangeTensorElement<5>;
}
