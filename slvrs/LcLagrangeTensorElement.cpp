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

// construct sequencer for use in mapping functions
    int lower[NDIM], upper[NDIM];
    for (unsigned d=0; d<NDIM; ++d) 
    {
      lower[d] = 0;
      upper[d] = numNodes[d];
    }
    nodeSeq = Lucee::RowMajorSequencer<NDIM>(Lucee::Region<NDIM, int>(lower, upper));

// initialize number of local nodes
    unsigned nlocal = basisCalc.getNumNodes();
    this->setNumNodes(nlocal);

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
      local2Global = Lucee::RowMajorIndexer<NDIM>(
        Lucee::Region<NDIM, int>(lower, upper));
    }
    else
    {
      int lower[NDIM], upper[NDIM];
      for (unsigned d=0; d<NDIM; ++d)
      {
        lower[d] = 0;
        upper[d] = (numNodes[d]-1)*gridRgn.getShape(d) + 1;
      }
      local2Global = Lucee::RowMajorIndexer<NDIM>(
        Lucee::Region<NDIM, int>(lower, upper));
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

    double dx[NDIM], vol2 = 1.0;
    for (unsigned d=0; d<NDIM; ++d)
    {
      dx[d] = grid.getDx(d);
      vol2 *= 0.5*dx[d];
    }

// mass matrix
    mass = Matrix<double>(nlocal, nlocal);
    basisCalc.getMassMatrix(mass);
    mass *= vol2;

// grad-stiff matrices
    for (unsigned d=0; d<NDIM; ++d)
    {
      gradStiff[d] = Matrix<double>(nlocal, nlocal);
      basisCalc.getGradStiffnessMatrix(d, gradStiff[d]);
      gradStiff[d] *= 2.0/dx[d]*vol2;
    }

// TODO: STIFFNESS MATRIX

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
  LagrangeTensorElement<NDIM>::getExclusiveNodeIndices(std::vector<unsigned>& ndIds)
  {
    basisCalc.getExclusiveNodeIndices(ndIds);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getLocalToGlobal(std::vector<int>& lgMap) const
  {
    int idx[NDIM];
    unsigned count = 0;
    nodeSeq.reset();
    while (nodeSeq.step())
    {
      nodeSeq.fillWithIndex(idx);
      for (unsigned d=0; d<NDIM; ++d)
        idx[d] += this->currIdx[d];
      lgMap[count++] = local2Global.getIndex(idx);
    }
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getSurfLowerLocalToGlobal(unsigned dir, std::vector<int>& lgMap) const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getSurfLowerLocalToGlobal(dir, lgMap);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getSurfUpperLocalToGlobal(unsigned dir, std::vector<int>& lgMap) const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getSurfUpperLocalToGlobal(dir, lgMap);
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
    return Lucee::NodalFiniteElementIfc<NDIM>::getWeights(w);
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
    return Lucee::NodalFiniteElementIfc<NDIM>::getLowerFaceMassMatrix(dir, NjNk);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getUpperFaceMassMatrix(unsigned dir, Lucee::Matrix<double>& NjNk) const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getUpperFaceMassMatrix(dir, NjNk);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getStiffnessMatrix(Lucee::Matrix<double>& DNjDNk) const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getStiffnessMatrix(DNjDNk);
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
    return Lucee::NodalFiniteElementIfc<NDIM>::getNumGaussNodes();
  }

  template <unsigned NDIM>
  unsigned
  LagrangeTensorElement<NDIM>::getNumSurfGaussNodes() const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getNumSurfGaussNodes();
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getGaussQuadData(Lucee::Matrix<double>& interpMat,
    Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getGaussQuadData(interpMat, ordinates, weights);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getSurfLowerGaussQuadData(unsigned dir, Lucee::Matrix<double>& interpMat,
    Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getSurfLowerGaussQuadData(dir, interpMat, ordinates, weights);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getSurfUpperGaussQuadData(unsigned dir, Lucee::Matrix<double>& interpMat,
    Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>:: getSurfUpperGaussQuadData(dir, interpMat, ordinates, weights);
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
    return Lucee::NodalFiniteElementIfc<NDIM>::extractFromField(fld, data);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::copyAllDataFromField(const Lucee::Field<NDIM, double>& fld, double *data)
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::copyAllDataFromField(fld, data);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::copyAllDataToField(const double *data, Lucee::Field<NDIM, double>& fld)
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::copyAllDataToField(data, fld);
  }

// instantiations
  template class LagrangeTensorElement<1>;
  template class LagrangeTensorElement<2>;
  template class LagrangeTensorElement<3>;
  template class LagrangeTensorElement<4>;
  template class LagrangeTensorElement<5>;
}
