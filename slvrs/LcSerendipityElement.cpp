/**
 * @file	LcSerendipityElement.cpp
 *
 * @brief Serendipity element implemented so far for 2 and 3 dimensions.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcExcept.h>
#include <LcGridIfc.h>
#include <LcSerendipityElement.h>

namespace Lucee
{
  using namespace Eigen;

// set module name
  template <> const char *SerendipityElement<1>::id = "SerendipityElement";
  template <> const char *SerendipityElement<2>::id = "SerendipityElement";
  template <> const char *SerendipityElement<3>::id = "SerendipityElement";
  template <> const char *SerendipityElement<4>::id = "SerendipityElement";
  template <> const char *SerendipityElement<5>::id = "SerendipityElement";

  template <unsigned NDIM>
  SerendipityElement<NDIM>::SerendipityElement()
    : NodalFiniteElementIfc<NDIM>(1)
  {

  }
  
  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::NodalFiniteElementIfc<NDIM>::readInput(tbl);
    polyOrder = (unsigned) tbl.getNumber("polyOrder");
    // Compute maximum polynomial power we want (no transformations)
    // This will be the maximum power of a single variable in anything
    // evaluated (e.g. product of basis functions)
    // Note: don't change this without paying attention to how maxPower is
    // used in the rest of the code! Affects more than just gaussPoints.
    if (polyOrder == 1)
      maxPower = 4;
    else
      maxPower = 3*polyOrder;
    
    if (NDIM == 1)
    {
      if (polyOrder < 4)
        this->setNumNodes(getSerendipityDimension(polyOrder, NDIM));
      else
      {
        Lucee::Except lce("SerendipityElement: Degree must be 1, 2, or 3.");
        lce << " Provided " << polyOrder << " instead";
        throw lce;
      }
      setupMatrices();
    }
    else if (NDIM == 2)
    {
      if (polyOrder < 4)
        this->setNumNodes(getSerendipityDimension(polyOrder, NDIM));
      else
      {
        Lucee::Except lce("SerendipityElement: Degree must be 1, 2, or 3.");
        lce << " Provided " << polyOrder << " instead";
        throw lce;
      }
      setupMatrices();
    }
    else if (NDIM == 3 || NDIM == 4 || NDIM == 5)
    {
      if (polyOrder < 3)
        this->setNumNodes(getSerendipityDimension(polyOrder, NDIM));
      else
      {
        Lucee::Except lce("SerendipityElement: Degree must be 1 or 2");
        lce << " Provided " << polyOrder << " instead";
        throw lce;
      }
      setupMatrices();
    }
    else
    {
      Lucee::Except lce("SerendipityElement: NDIM not supported.");
      lce << " Provided " << NDIM << " instead";
      throw lce;
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getExclusiveNodeIndices(std::vector<int>& ndIds)
  {
    if (NDIM == 1)
    {
      ndIds.clear();
      ndIds.resize(polyOrder);

      for (int nodeIndex = 0; nodeIndex < polyOrder; nodeIndex++)
        ndIds[nodeIndex] = nodeIndex;
    }
    else if (NDIM == 2)
    {
      ndIds.clear();
      if (polyOrder == 1)
      {
        ndIds.resize(1);
        ndIds[0] = 0;
      }
      else if (polyOrder == 2)
      {
        ndIds.resize(3);
        ndIds[0] = 0;
        ndIds[1] = 1;
        ndIds[2] = 7;
      }
      else if (polyOrder == 3)
      {
        ndIds.resize(5);
        ndIds[0] = 0;
        ndIds[1] = 1;
        ndIds[2] = 2;
        ndIds[3] = 11;
        ndIds[4] = 10;
      }
      else if (polyOrder == 4)
      {
        ndIds.resize(8);
        ndIds[0] = 0;
        ndIds[1] = 1;
        ndIds[2] = 2;
        ndIds[3] = 3;
        ndIds[4] = 15;
        ndIds[5] = 14;
        ndIds[6] = 13;
        ndIds[7] = 16;
      }
    }
    else
    {
      Lucee::Except lce("SerendipityElement::getExclusiveNodeIndices NDIM must be 1 or 2.");
      lce << " Provided " << NDIM << " instead";
      throw lce;
    }
  }

  template <unsigned NDIM>
  unsigned
  SerendipityElement<NDIM>::getNumSurfLowerNodes(unsigned dir) const
  {
    return getSerendipityDimension(polyOrder, NDIM-1);
  }

  template <unsigned NDIM>
  unsigned
  SerendipityElement<NDIM>::getNumSurfUpperNodes(unsigned dir) const
  {
    return getSerendipityDimension(polyOrder, NDIM-1);
  }

  template <unsigned NDIM>
  unsigned
  SerendipityElement<NDIM>::getNumGlobalNodes() const
  {
    if (NDIM != 1)
      throw Lucee::Except("SerendipityElement::getNumGlobalNodes: Not implemented for NDIM > 1!");
    // get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->template getGrid<Lucee::StructuredGridBase<NDIM> >();
    Lucee::Region<NDIM, int> gridRgn = grid.getGlobalRegion();

    int numGlobalNodes = 1 + polyOrder*gridRgn.getShape(0);
    return numGlobalNodes;
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getLocalToGlobal(std::vector<int>& lgMap) const
  {
    std::cout << "this->currIdx[0] = " << this->currIdx[0] << std::endl;

    if (NDIM == 1)
    {
    // Loop over excusive nodes
    for (int nodeIndex = 0; nodeIndex < polyOrder; nodeIndex++)
      lgMap[nodeIndex] = this->currIdx[0]*polyOrder + nodeIndex;
    }
    else throw Lucee::Except("SerendipityElement::getLocalToGlobal: Not implemented!");
    /*
    // TODO: untested, not sure what this does anymore.
    // determine number of global degrees of freedom
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->template getGrid<Lucee::StructuredGridBase<NDIM> >();
    Lucee::Region<NDIM, int> grgn = grid.getGlobalRegion();
    unsigned cellsPerDim[NDIM];
    unsigned currIndex[NDIM];

    for (unsigned dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      cellsPerDim[dimIndex] = grgn.getShape(dimIndex);
      // Get global index of current cell, too
      currIndex[dimIndex]   = this->currIdx[dimIndex];
    }

    // Number of nodes on an edge
    unsigned nodesPerSide = polyOrder + 1;
    // Number of nodes on a 2-D serendipity element of same degree
    unsigned nodesPerCell2D;
    // Figure out number of nodes per row
    std::vector<unsigned> nodesPerRow2D(nodesPerSide);
    std::vector<unsigned> globalNodesPerRow2D(nodesPerSide);
    if (polyOrder == 1)
    {
      nodesPerCell2D = 4;
      nodesPerRow2D[0] = 2;
      nodesPerRow2D[1] = 2;
    }
    else if (polyOrder == 2)
    {
      nodesPerCell2D = 8;
      nodesPerRow2D[0] = 3;
      nodesPerRow2D[1] = 2;
      nodesPerRow2D[2] = 3;
    }
    else if (polyOrder == 3)
    {
      nodesPerCell2D = 12;
      nodesPerRow2D[0] = 4;
      nodesPerRow2D[1] = 2;
      nodesPerRow2D[2] = 2;
      nodesPerRow2D[3] = 4;
    }

    // Total number of nodes in one cell row of the mesh, minus the 'top row'
    // Type A nodes are those on the top and bottom rows of an element
    unsigned typeAPad = 0;
    for (unsigned rowIndex = 0; rowIndex < nodesPerSide; rowIndex++)
    {
      globalNodesPerRow2D[rowIndex] = (nodesPerRow2D[rowIndex]-1)*cellsPerDim[0] + 1;
      if (rowIndex != nodesPerSide-1)
        typeAPad += globalNodesPerRow2D[rowIndex];
    }

    std::vector<unsigned> nodeIndices(nodesPerCell2D);*/
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getSurfLowerLocalToGlobal(unsigned dim,
    std::vector<int>& lgMap) const
  {
    if (NDIM == 1)
      lgMap[0] = this->currIdx[0]*polyOrder;
    else throw Lucee::Except("SerendipityElement::getSurfLowerLocalToGlobal: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getSurfUpperLocalToGlobal(unsigned dim,
    std::vector<int>& lgMap) const
  {
    if (NDIM == 1)
      lgMap[0] = (this->currIdx[0]+1)*polyOrder;
    else throw Lucee::Except("SerendipityElement::getSurfUpperLocalToGlobal: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getSurfLowerNodeNums(unsigned dir,
    std::vector<int>& nodeNum) const
  { 
    if (NDIM == 1)
      nodeNum[0] = 0;
    else
    {
      int counter = 0;

      for (int i = 0; i < nodeList.rows(); i++)
      {
        if (nodeList(i,dir) == -1)
        {
          nodeNum[counter] = i;
          counter++;
        }
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getSurfUpperNodeNums(unsigned dir,
    std::vector<int>& nodeNum) const
  {
    if (NDIM == 1)
      nodeNum[0] = polyOrder;
    else
    {
      int counter = 0;

      for (int i = 0; i < nodeList.rows(); i++)
      {
        if (nodeList(i,dir) == 1)
        {
          nodeNum[counter] = i;
          counter++;
        }
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getNodalCoordinates(Lucee::Matrix<double>& nodeCoords)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->template getGrid<Lucee::StructuredGridBase<NDIM> >();
    // Set index and get centroid coordinate
    grid.setIndex(this->currIdx);
    double xc[NC];

    grid.getCentroid(xc);

    // Loop over all node locations on reference element and convert them
    // to appropriate coordinates
    for (int i = 0; i < this->getNumNodes(); i++)
    {
      for (int dim = 0; dim < NC; dim++)
      {
        if (dim < NDIM)
          nodeCoords(i, dim) = xc[dim] + nodeList(i, dim)*0.5*dq[dim];
        else
          nodeCoords(i, dim) = 0;
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getWeights(std::vector<double>& w)
  {
    // Todo (2)
    throw Lucee::Except("SerendipityElement::getWeights: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getSurfUpperWeights(unsigned dir, std::vector<double>& w)
  {
    // Todo (2)
    throw Lucee::Except("SerendipityElement::getSurfUpperWeights: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getSurfLowerWeights(unsigned dir, std::vector<double>& w)
  {
    // Todo (2)
    throw Lucee::Except("SerendipityElement::getSurfLowerWeights: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getMassMatrix(Lucee::Matrix<double>& NjNk) const
  {
    for (int i = 0; i < refMass.rows(); i++)
      for (int j = 0; j < refMass.cols(); j++)
        NjNk(i, j) = refMass(i, j);
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getLowerFaceMassMatrix(unsigned dir,
    Lucee::Matrix<double>& NjNk) const
  {
    for (int i = 0; i < refFaceMassLower[dir].rows(); i++)
      for (int j = 0; j < refFaceMassLower[dir].cols(); j++)
        NjNk(i, j) = refFaceMassLower[dir](i, j);
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getUpperFaceMassMatrix(unsigned dir,
    Lucee::Matrix<double>& NjNk) const
  {
    for (int i = 0; i < refFaceMassUpper[dir].rows(); i++)
      for (int j = 0; j < refFaceMassUpper[dir].cols(); j++)
        NjNk(i, j) = refFaceMassUpper[dir](i, j);
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getStiffnessMatrix(Lucee::Matrix<double>& DNjDNk) const
  {
    for (int i = 0; i < refStiffness.rows(); i++)
      for (int j = 0; j < refStiffness.cols(); j++)
        DNjDNk(i, j) = refStiffness(i, j);
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getGradStiffnessMatrix(
    unsigned dir, Lucee::Matrix<double>& DNjNk) const
  {
    for (int i = 0; i < refGradStiffness[dir].rows(); i++)
      for (int j = 0; j < refGradStiffness[dir].cols(); j++)
        DNjNk(i, j) = refGradStiffness[dir](i, j);
  }

  template <unsigned NDIM>
  unsigned
  SerendipityElement<NDIM>::getNumGaussNodes() const
  {
    return gaussNodeList.rows();
  }

  template <unsigned NDIM>
  unsigned
  SerendipityElement<NDIM>::getNumSurfGaussNodes() const
  {
    return gaussNodeListLowerSurf[0].rows();
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getGaussQuadData(Lucee::Matrix<double>& interpMat,
    Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const
  {
    double weightScale = 1.0;

    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
      weightScale *= 0.5*dq[dimIndex];

    for (int gaussIndex = 0; gaussIndex < gaussNodeList.rows(); gaussIndex++)
    {
      if (NDIM == 1)
      {
        ordinates(gaussIndex,0) = gaussNodeList(gaussIndex,0);
        ordinates(gaussIndex,1) = 0;
        ordinates(gaussIndex,2) = 0;
      }
      else if (NDIM == 2)
      {
        ordinates(gaussIndex,0) = gaussNodeList(gaussIndex,0);
        ordinates(gaussIndex,1) = gaussNodeList(gaussIndex,1);
        ordinates(gaussIndex,2) = 0;
      }
      else
      {
        for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
          ordinates(gaussIndex,dimIndex) = gaussNodeList(gaussIndex,dimIndex);
      }

      weights[gaussIndex] = gaussNodeList(gaussIndex,NDIM)*weightScale;

      // Copy interpolation matrix
      for (int basisIndex = 0; basisIndex < functionEvaluations.cols(); basisIndex++)
        interpMat(gaussIndex, basisIndex) = functionEvaluations(gaussIndex, basisIndex);
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getSurfLowerGaussQuadData(unsigned dir, Lucee::Matrix<double>& interpMat,
    Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const
  {
    double weightScale = 1.0;

    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      if (dimIndex != dir)
        weightScale *= 0.5*dq[dimIndex];
    }

    for (int gaussIndex = 0; gaussIndex < gaussNodeListLowerSurf[dir].rows(); gaussIndex++)
    {
      if (NDIM == 1)
      {
        ordinates(gaussIndex,0) = gaussNodeListLowerSurf[dir](gaussIndex,0);
        ordinates(gaussIndex,1) = 0;
        ordinates(gaussIndex,2) = 0;
      }
      else if (NDIM == 2)
      {
        ordinates(gaussIndex,0) = gaussNodeListLowerSurf[dir](gaussIndex,0);
        ordinates(gaussIndex,1) = gaussNodeListLowerSurf[dir](gaussIndex,1);
        ordinates(gaussIndex,2) = 0;
      }
      else
      {
        for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
          ordinates(gaussIndex,dimIndex) = gaussNodeListLowerSurf[dir](gaussIndex,dimIndex);
      }

      weights[gaussIndex] = gaussNodeListLowerSurf[dir](gaussIndex,NDIM)*weightScale;

      // Copy interpolation matrix
      for (int basisIndex = 0; basisIndex < lowerSurfaceEvaluations[dir].cols(); basisIndex++)
        interpMat(gaussIndex, basisIndex) = lowerSurfaceEvaluations[dir](gaussIndex, basisIndex);
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getSurfUpperGaussQuadData(unsigned dir, Lucee::Matrix<double>& interpMat,
    Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const
  {
    double weightScale = 1.0;

    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      if (dimIndex != dir)
        weightScale *= 0.5*dq[dimIndex];
    }

    for (int gaussIndex = 0; gaussIndex < gaussNodeListUpperSurf[dir].rows(); gaussIndex++)
    {
      if (NDIM == 1)
      {
        ordinates(gaussIndex,0) = gaussNodeListUpperSurf[dir](gaussIndex,0);
        ordinates(gaussIndex,1) = 0;
        ordinates(gaussIndex,2) = 0;
      }
      else if (NDIM == 2)
      {
        ordinates(gaussIndex,0) = gaussNodeListUpperSurf[dir](gaussIndex,0);
        ordinates(gaussIndex,1) = gaussNodeListUpperSurf[dir](gaussIndex,1);
        ordinates(gaussIndex,2) = 0;
      }
      else
      {
        for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
          ordinates(gaussIndex,dimIndex) = gaussNodeListUpperSurf[dir](gaussIndex,dimIndex);
      }

      weights[gaussIndex] = gaussNodeListUpperSurf[dir](gaussIndex,NDIM)*weightScale;

      // Copy interpolation matrix
      for (int basisIndex = 0; basisIndex < upperSurfaceEvaluations[dir].cols(); basisIndex++)
        interpMat(gaussIndex, basisIndex) = upperSurfaceEvaluations[dir](gaussIndex, basisIndex);
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getMomentMatrix(unsigned p, Lucee::Matrix<double>& momMat) const
  {
    if (p > maxMoment || NDIM != 2)
    {
      // moments higher than 3 not supported for now
      Lucee::Except lce("SerendipityElement::getMomentMatrix: Moment matrix of order ");
      lce << p << " not supported. NDIM must also be 2.";
      throw lce;
    }
    else
    {
      for (int rowIndex = 0; rowIndex < momMatrix[p].rows(); rowIndex++)
        for (int colIndex = 0; colIndex < momMatrix[p].cols(); colIndex++)
          momMat(rowIndex,colIndex) = momMatrix[p](rowIndex,colIndex);
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getDiffusionMatrices(std::vector<Lucee::Matrix<double> >& iMat_o,
    std::vector<Lucee::Matrix<double> >& lowerMat_o, std::vector<Lucee::Matrix<double> >& upperMat_o) const
  {
    if (polyOrder > 2)
      Lucee::Except("SerendipityElement::getDiffusionMatrices: Not implemented for higher than quadratic!");
    if (NDIM != 2)
      Lucee::Except("SerendipityElement::getDiffusionMatrices: Only implemented for 2D");

    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      for (int rowIndex = 0; rowIndex < iMatDiffusion[dimIndex].rows(); rowIndex++)
      {
        for (int colIndex = 0; colIndex < iMatDiffusion[dimIndex].cols(); colIndex++)
        {
          iMat_o[dimIndex](rowIndex,colIndex) = iMatDiffusion[dimIndex](rowIndex,colIndex);
          lowerMat_o[dimIndex](rowIndex,colIndex) = lowerMatDiffusion[dimIndex](rowIndex,colIndex);
          upperMat_o[dimIndex](rowIndex,colIndex) = upperMatDiffusion[dimIndex](rowIndex,colIndex);
        }
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getHyperDiffusionMatrices(std::vector<Lucee::Matrix<double> >& iMat_o,
    std::vector<Lucee::Matrix<double> >& lowerMat_o, std::vector<Lucee::Matrix<double> >& upperMat_o) const
  {
    if (polyOrder > 2)
      Lucee::Except("SerendipityElement::getHyperDiffusionMatrices: Not implemented for higher than quadratic!");
    if (NDIM != 2)
      Lucee::Except("SerendipityElement::getHyperDiffusionMatrices: Only implemented for 2D");

    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      for (int rowIndex = 0; rowIndex < iMatHyperDiffusion[dimIndex].rows(); rowIndex++)
      {
        for (int colIndex = 0; colIndex < iMatHyperDiffusion[dimIndex].cols(); colIndex++)
        {
          iMat_o[dimIndex](rowIndex,colIndex) = iMatHyperDiffusion[dimIndex](rowIndex,colIndex);
          lowerMat_o[dimIndex](rowIndex,colIndex) = lowerMatHyperDiffusion[dimIndex](rowIndex,colIndex);
          upperMat_o[dimIndex](rowIndex,colIndex) = upperMatHyperDiffusion[dimIndex](rowIndex,colIndex);
        }
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getLowerReflectingBcMapping(unsigned dir,
        std::vector<unsigned>& nodeMap) const
  {
    if (polyOrder > 2)
      Lucee::Except("SerendipityElement::getLowerReflectingBcMapping: Not implemented for higher than quadratic!");
    if (NDIM != 2)
      Lucee::Except("SerendipityElement::getLowerReflectingBcMapping: Only implemented for 2D");

    if(dir == 0)
    {
      if (polyOrder == 1)
      {
        nodeMap[0] = 1;
        nodeMap[1] = 0;
        nodeMap[2] = 3;
        nodeMap[3] = 2;
      }
      else if (polyOrder == 2)
      {
        nodeMap[0] = 2;
        nodeMap[1] = 1;
        nodeMap[2] = 0;
        nodeMap[3] = 7;
        nodeMap[4] = 6;
        nodeMap[5] = 5;
        nodeMap[6] = 4;
        nodeMap[7] = 3;
      }
    }
    else if (dir == 1)
    {
      if (polyOrder == 1)
      {
        nodeMap[0] = 3;
        nodeMap[1] = 2;
        nodeMap[2] = 1;
        nodeMap[3] = 0;
      }
      else if (polyOrder == 2)
      {
        nodeMap[0] = 6;
        nodeMap[1] = 5;
        nodeMap[2] = 4;
        nodeMap[3] = 3;
        nodeMap[4] = 2;
        nodeMap[5] = 1;
        nodeMap[6] = 0;
        nodeMap[7] = 7;
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getUpperReflectingBcMapping(unsigned dir,
        std::vector<unsigned>& nodeMap) const
  {
    if (polyOrder > 2)
      Lucee::Except("SerendipityElement::getUpperReflectingBcMapping: Not implemented for higher than quadratic!");
    if (NDIM != 2)
      Lucee::Except("SerendipityElement::getUpperReflectingBcMapping: Only implemented for 2D");

    if(dir == 0)
    {
      if (polyOrder == 1)
      {
        nodeMap[0] = 1;
        nodeMap[1] = 0;
        nodeMap[2] = 3;
        nodeMap[3] = 2;
      }
      else if (polyOrder == 2)
      {
        nodeMap[0] = 2;
        nodeMap[1] = 1;
        nodeMap[2] = 0;
        nodeMap[3] = 7;
        nodeMap[4] = 6;
        nodeMap[5] = 5;
        nodeMap[6] = 4;
        nodeMap[7] = 3;
      }
    }
    else if (dir == 1)
    {
      if (polyOrder == 1)
      {
        nodeMap[0] = 3;
        nodeMap[1] = 2;
        nodeMap[2] = 1;
        nodeMap[3] = 0;
      }
      else if (polyOrder == 2)
      {
        nodeMap[0] = 6;
        nodeMap[1] = 5;
        nodeMap[2] = 4;
        nodeMap[3] = 3;
        nodeMap[4] = 2;
        nodeMap[5] = 1;
        nodeMap[6] = 0;
        nodeMap[7] = 7;
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getLowerFaceToInteriorMapping(unsigned dir,
        Lucee::Matrix<double>& faceToIntMap) const
  {
    if (polyOrder > 3)
      Lucee::Except("SerendipityElement::getLowerFaceToInteriorMapping: Not implemented for higher than cubic!");
    if (NDIM != 2)
      Lucee::Except("SerendipityElement::getLowerFaceToInteriorMapping: Only implemented for 2D");

    for (int rowIndex = 0; rowIndex < lowerFaceToInteriorMapMatrices[dir].rows(); rowIndex++)
      for (int colIndex = 0; colIndex < lowerFaceToInteriorMapMatrices[dir].cols(); colIndex++)
        faceToIntMap(rowIndex, colIndex) = lowerFaceToInteriorMapMatrices[dir](rowIndex, colIndex);
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::extractFromField(const Lucee::Field<NDIM, double>& fld,
    std::vector<double>& data)
  {
    // Only implemented for 2-D elements
    if (NDIM != 2)
    {
      Lucee::Except lce("SerendipityElement::extractFromField: NDIM ");
      lce << NDIM << " not supported";
      throw lce;
    }

    Lucee::ConstFieldPtr<double> fldPtr = fld.createConstPtr();
    Lucee::ConstFieldPtr<double> fldPtr_r = fld.createConstPtr();
    Lucee::ConstFieldPtr<double> fldPtr_t = fld.createConstPtr();
    Lucee::ConstFieldPtr<double> fldPtr_tr = fld.createConstPtr();

    // attach pointers to proper locations
    fld.setPtr(fldPtr, this->currIdx[0], this->currIdx[1]);
    fld.setPtr(fldPtr_r, this->currIdx[0]+1, this->currIdx[1]);
    fld.setPtr(fldPtr_t, this->currIdx[0], this->currIdx[1]+1);
    fld.setPtr(fldPtr_tr, this->currIdx[0]+1, this->currIdx[1]+1);

    // copy data
    if (polyOrder == 1)
    {
      data[0] = fldPtr[0];
      data[1] = fldPtr_r[0];
      data[2] = fldPtr_tr[0];
      data[3] = fldPtr_t[0];
    }
    else if (polyOrder == 2)
    {
      data[0] = fldPtr[0];
      data[1] = fldPtr[1];
      data[2] = fldPtr_r[0];
      data[3] = fldPtr_r[2];
      data[4] = fldPtr_tr[0];
      data[5] = fldPtr_t[1];
      data[6] = fldPtr_t[0];
      data[7] = fldPtr[2];
    }
    else if (polyOrder == 3)
    {
      data[0] = fldPtr[0];
      data[1] = fldPtr[1];
      data[2] = fldPtr[2];
      data[3] = fldPtr_r[0];
      data[4] = fldPtr_r[3];
      data[5] = fldPtr_r[4];
      data[6] = fldPtr_tr[0];
      data[7] = fldPtr_t[2];
      data[8] = fldPtr_t[1];
      data[9] = fldPtr_t[0];
      data[10] = fldPtr[4];
      data[11] = fldPtr[3];
    }
    else if (polyOrder == 4)
    {
      data[0] = fldPtr[0];
      data[1] = fldPtr[1];
      data[2] = fldPtr[2];
      data[3] = fldPtr[3];
      data[4] = fldPtr_r[0];
      data[5] = fldPtr_r[4];
      data[6] = fldPtr_r[5];
      data[7] = fldPtr_r[6];
      data[8] = fldPtr_tr[0];
      data[9] = fldPtr_t[3];
      data[10] = fldPtr_t[2];
      data[11] = fldPtr_t[1];
      data[12] = fldPtr_t[0];
      data[13] = fldPtr[6];
      data[14] = fldPtr[5];
      data[15] = fldPtr[4];
      data[16] = fldPtr[7];
    }
    else
    {
      Lucee::Except lce("SerendipityElement::extractFromField: polyOrder ");
      lce << polyOrder << " not supported";
      throw lce;
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::copyAllDataFromField(const Lucee::Field<NDIM, double>& fld,
    double *data)
  {
    throw Lucee::Except("SerendipityElement::copyAllDataFromField: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::copyAllDataToField(const double *data, 
    Lucee::Field<NDIM, double>& fld)
  {
    throw Lucee::Except("SerendipityElement::copyAllDataToField: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::evalBasis(double xc[NDIM], std::vector<double>& vals) const
  {
    // Create vector from xc
    Eigen::VectorXd nodeVec(NDIM);
    for (int componentIndex = 0; componentIndex < NDIM; componentIndex++)
      nodeVec(componentIndex) = xc[componentIndex];

    for (int basisIndex = 0; basisIndex < functionVector.size(); basisIndex++)
      vals[basisIndex] = evalPolynomial(functionVector[basisIndex], nodeVec);
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::setupMatrices()
  {
    // Get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->template getGrid<Lucee::StructuredGridBase<NDIM> >();
    // Get grid spacing
    for (int i = 0; i < NDIM; i++)
    {
      dq[i] = grid.getDx(i);
      dq2[i] = dq[i]*dq[i];
    }

    // Populate nodeList according to polyOrder
    nodeList = Eigen::MatrixXd(this->getNumNodes(),NDIM);
    getNodeList(nodeList, polyOrder, NDIM);

    //std::vector<blitz::Array<double,NDIM> > functionVector;
    computeBasisFunctions(functionVector, nodeList, polyOrder);

    // Compute gaussian quadrature weights and locations in 1-D
    numGaussPoints = (unsigned)((maxPower+1)/2.0 + 0.5);
    gaussPoints  = std::vector<double>(numGaussPoints);
    gaussWeights = std::vector<double>(numGaussPoints);
    legendre_set(numGaussPoints, &gaussPoints[0], &gaussWeights[0]);
    
    // Figure out how many gaussian points there are
    int totalVolumeGaussNodes = 1;

    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
      totalVolumeGaussNodes = gaussPoints.size()*totalVolumeGaussNodes;

    int totalSurfaceGaussNodes = totalVolumeGaussNodes/gaussPoints.size();
    
    gaussNodeList = Eigen::MatrixXd::Zero(totalVolumeGaussNodes, NDIM+1);
    functionEvaluations  = Eigen::MatrixXd(totalVolumeGaussNodes, this->getNumNodes());
    // 3D Array containing derivatives of basis functions (cols) evaluated at gaussian
    // integration locations (rows). Correspondance between column and gaussian node set
    // is kept track of in gaussNodeList
    blitz::Array<double,3> functionDEvaluations(totalVolumeGaussNodes, this->getNumNodes(),NDIM);
    
    gaussNodeListUpperSurf  = std::vector<Eigen::MatrixXd>(NDIM);
    gaussNodeListLowerSurf  = std::vector<Eigen::MatrixXd>(NDIM);
    upperSurfaceEvaluations = std::vector<Eigen::MatrixXd>(NDIM);
    lowerSurfaceEvaluations = std::vector<Eigen::MatrixXd>(NDIM);

    Eigen::MatrixXd gaussNodeListSurf = Eigen::MatrixXd::Zero(totalSurfaceGaussNodes, NDIM);

    if (NDIM > 1)
    {
      // TODO: roll surface and volume quadrature code into one function
      // since they do the same thing, but in diff dimensions
      // Compute surface gaussian quadrature locations
      int surfShape[NDIM-1];
      int surfIdx[NDIM-1];
      
      for (int dimIndex = 0; dimIndex < NDIM-1; dimIndex++)
        surfShape[dimIndex] = gaussPoints.size();

      Lucee::Region<NDIM-1, int> surfRegion(surfShape);

      Lucee::RowMajorSequencer<NDIM-1> surfSeq = RowMajorSequencer<NDIM-1>(surfRegion);
      Lucee::RowMajorIndexer<NDIM-1> surfIdxr = RowMajorIndexer<NDIM-1>(surfRegion);

      
      // Find all quadrature locations on a NDIM-1 surface
      while(surfSeq.step())
      {
        surfSeq.fillWithIndex(surfIdx);
        int nodeNumber = surfIdxr.getIndex(surfIdx);

        gaussNodeListSurf(nodeNumber, NDIM-1) = 1.0;
        for (int dimIndex = 0; dimIndex < NDIM-1; dimIndex++)
        {
          gaussNodeListSurf(nodeNumber, dimIndex) = gaussPoints[surfIdx[dimIndex]];
          gaussNodeListSurf(nodeNumber, NDIM-1)    *= gaussWeights[surfIdx[dimIndex]];
        }
      }
    }
    else
    {
      gaussNodeListSurf(0, 0) = 1.0;
    }

    // Evaluate quadrature points on all surfaces
    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      // Initialize matrices
      gaussNodeListUpperSurf[dimIndex]  = Eigen::MatrixXd::Zero(totalSurfaceGaussNodes, NDIM+1);
      gaussNodeListLowerSurf[dimIndex]  = Eigen::MatrixXd::Zero(totalSurfaceGaussNodes, NDIM+1);
      upperSurfaceEvaluations[dimIndex] = Eigen::MatrixXd::Zero(totalSurfaceGaussNodes, functionVector.size());
      lowerSurfaceEvaluations[dimIndex] = Eigen::MatrixXd::Zero(totalSurfaceGaussNodes, functionVector.size());
      
      Eigen::VectorXd gaussNodeVec = Eigen::VectorXd::Zero(NDIM);
      // Evaluate all basis functions at all upper and lower surfaces
      for (int nodeIndex = 0; nodeIndex < gaussNodeListSurf.rows(); nodeIndex++)
      {
        int rollingIndex = 0;
        
        for (int coordIndex = 0; coordIndex < NDIM; coordIndex++)
        {
          if (coordIndex == dimIndex)
            gaussNodeVec(dimIndex) = 1;
          else
          {
            gaussNodeVec(coordIndex) = gaussNodeListSurf(nodeIndex, rollingIndex);
            rollingIndex++;
          }
          gaussNodeListUpperSurf[dimIndex](nodeIndex, coordIndex) = gaussNodeVec(coordIndex);
        }
        gaussNodeListUpperSurf[dimIndex](nodeIndex, NDIM) = gaussNodeListSurf(nodeIndex, NDIM-1);

        // Evaluate all basis functions at this node location
        for (int basisIndex = 0; basisIndex < functionVector.size(); basisIndex++)
          upperSurfaceEvaluations[dimIndex](nodeIndex, basisIndex) = evalPolynomial(functionVector[basisIndex],gaussNodeVec);
        
        // Copy upper surface matrix into lower one first, then replace surface index
        gaussNodeListLowerSurf[dimIndex].row(nodeIndex) = gaussNodeListUpperSurf[dimIndex].row(nodeIndex);
        gaussNodeListLowerSurf[dimIndex](nodeIndex, dimIndex) = -1;
        gaussNodeVec(dimIndex) = -1;
        // Evaluate all basis functions at this node location
        for (int basisIndex = 0; basisIndex < functionVector.size(); basisIndex++)
          lowerSurfaceEvaluations[dimIndex](nodeIndex, basisIndex) = evalPolynomial(functionVector[basisIndex],gaussNodeVec);
      }
    }

    // Evaluate all basis functions and their derivatives in every direction at all gaussian volume nodes
    // First compute all volume gaussian quadrature locations
    int volShape[NDIM];
    int volIdx[NDIM];
    
    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
      volShape[dimIndex] = gaussPoints.size();

    Lucee::Region<NDIM, int> volRegion(volShape);

    // Needs col-major right now to help project 1-D data into 2-D element
    Lucee::ColMajorSequencer<NDIM> volSeq = ColMajorSequencer<NDIM>(volRegion);
    Lucee::ColMajorIndexer<NDIM> volIdxr = ColMajorIndexer<NDIM>(volRegion);

    // Find all quadrature locations on a NDIM volume
    while(volSeq.step())
    {
      volSeq.fillWithIndex(volIdx);
      int nodeNumber = volIdxr.getIndex(volIdx);

      gaussNodeList(nodeNumber, NDIM) = 1.0;
      Eigen::VectorXd gaussNodeVec = Eigen::VectorXd::Zero(NDIM);
      for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
      {
        gaussNodeVec(dimIndex)              = gaussPoints[volIdx[dimIndex]];
        gaussNodeList(nodeNumber, dimIndex) = gaussPoints[volIdx[dimIndex]];
        gaussNodeList(nodeNumber, NDIM)    *= gaussWeights[volIdx[dimIndex]];
      }

      for (int basisIndex = 0; basisIndex < functionVector.size(); basisIndex++)
      {
        functionEvaluations(nodeNumber, basisIndex) = evalPolynomial(functionVector[basisIndex],gaussNodeVec);
        
        for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
          functionDEvaluations(nodeNumber, basisIndex, dimIndex) = evalPolynomial(computePolynomialDerivative(functionVector[basisIndex],dimIndex),gaussNodeVec);
      }
    }

    // Resize+Initialize (most of) the output matrices we need
    resizeMatrices();

    // Call various functions to populate the matrices

    if (NDIM == 2)
    {
      // Comput moment matrices on reference element
      setupMomentMatrices();
      // Scale into physical space
      for(int pIndex = 0; pIndex < maxMoment+1; pIndex++)
        for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
          momMatrix[pIndex] *= 0.5*dq[dimIndex];
      
      // Set up diffusion matrices
      iMatDiffusion          = std::vector<Eigen::MatrixXd>(NDIM);
      iMatHyperDiffusion     = std::vector<Eigen::MatrixXd>(NDIM);
      lowerMatDiffusion      = std::vector<Eigen::MatrixXd>(NDIM);
      lowerMatHyperDiffusion = std::vector<Eigen::MatrixXd>(NDIM);
      upperMatDiffusion      = std::vector<Eigen::MatrixXd>(NDIM);
      upperMatHyperDiffusion = std::vector<Eigen::MatrixXd>(NDIM);

      for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
      {
        iMatDiffusion[dimIndex]          = Eigen::MatrixXd::Zero(functionVector.size(),functionVector.size());
        iMatHyperDiffusion[dimIndex]     = Eigen::MatrixXd::Zero(functionVector.size(),functionVector.size());
        lowerMatDiffusion[dimIndex]      = Eigen::MatrixXd::Zero(functionVector.size(),functionVector.size());
        lowerMatHyperDiffusion[dimIndex] = Eigen::MatrixXd::Zero(functionVector.size(),functionVector.size());
        upperMatDiffusion[dimIndex]      = Eigen::MatrixXd::Zero(functionVector.size(),functionVector.size());
        upperMatHyperDiffusion[dimIndex] = Eigen::MatrixXd::Zero(functionVector.size(),functionVector.size());
      }

      // Explicity assign values to matrix elements (very long)
      #include <LcSerendipityElementDiffusionOutput>
      #include <LcSerendipityElementHyperDiffusionOutput>

      // Set up 1D to 2D mapping matrix
      lowerFaceToInteriorMapMatrices = std::vector<Eigen::MatrixXd>(NDIM);
      
      for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
        lowerFaceToInteriorMapMatrices[dimIndex] = Eigen::MatrixXd::Zero(functionVector.size(), polyOrder + 1);

      // Explicitly assign values to matrix elements (somewhat long)
      #include <LcSerendipityElement2DFaceToInteriorOutput>
    }

    computeMass(refMass);
    computeStiffness(functionDEvaluations,refStiffness);

    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      computeFaceMass(dimIndex, refFaceMassLower[dimIndex], refFaceMassUpper[dimIndex]);
      computeGradStiffness(functionDEvaluations, dimIndex, refGradStiffness[dimIndex]);
    }

    // Scale the matrices computed on reference element into physical space
    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      refMass   *= 0.5*dq[dimIndex];
      refStiffness *= 0.5*dq[dimIndex];

      // Scale face-mass matrices
      for (int matrixIndex = 0; matrixIndex < NDIM; matrixIndex++)
      {
        if (dimIndex != matrixIndex)
        {
          refFaceMassUpper[matrixIndex] *= 0.5*dq[dimIndex];
          refFaceMassLower[matrixIndex] *= 0.5*dq[dimIndex];
          refGradStiffness[matrixIndex] *= 0.5*dq[dimIndex];
        }
      }
    }

    //computeTransformationScales();
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::resizeMatrices()
  {
    int numNodes = this->getNumNodes();

    refMass          = Eigen::MatrixXd::Zero(numNodes, numNodes);
    refStiffness     = Eigen::MatrixXd::Zero(numNodes, numNodes);

    refFaceMassLower = std::vector<Eigen::MatrixXd>(NDIM);
    refFaceMassUpper = std::vector<Eigen::MatrixXd>(NDIM);
    refGradStiffness = std::vector<Eigen::MatrixXd>(NDIM);

    for (unsigned dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      refFaceMassLower[dimIndex] = Eigen::MatrixXd::Zero(numNodes, getNumSurfLowerNodes(dimIndex));
      refFaceMassUpper[dimIndex] = Eigen::MatrixXd::Zero(numNodes, getNumSurfUpperNodes(dimIndex));
      refGradStiffness[dimIndex] = Eigen::MatrixXd::Zero(numNodes, numNodes);
    }
  }
 
  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getNodeList(Eigen::MatrixXd& nodeMatrix, int degree, int dimension)
  {
    if (dimension == 1)
    {
      if (degree == 1)
        nodeMatrix << -1,
                      1;
      else if (degree == 2)
        nodeMatrix << -1,
                       0,
                       1;
      else if (degree == 3)
        nodeMatrix << -1,
                      -1/3.0,
                      1/3.0,
                      1;
      else if (degree == 4)
        nodeMatrix << -1,
                      -0.5,
                      0.5,
                      1;
    }
    else if (dimension == 2)
    {
      if (degree == 1)
        nodeMatrix << -1,-1,
                      1,-1,
                      -1,1,
                      1,1;
      else if (degree == 2)
        nodeMatrix << -1,-1,
                    0,-1,
                    1,-1,
                    -1,0,
                    1,0,
                    -1,1,
                    0,1,
                    1,1;
      else if (degree == 3)
        nodeMatrix << -1,-1,
                    -1/3.0,-1,
                    1/3.0,-1,
                    1,-1,
                    -1,-1/3.0;
                    1,-1/3.0,
                    -1,1/3.0,
                    1,1/3.0,
                    -1,1,
                    -1/3.0,1,
                    1/3.0,1,
                    1,1;
    }
    else if (dimension == 3)
    {
      if (degree == 1)
        nodeMatrix << -1,-1,-1,
                     1,-1,-1,
                     -1,1,-1,
                     1,1,-1,
                     -1,-1,1,
                     1,-1,1,
                     -1,1,1,
                     1,1,1;
      else if (degree == 2)
        nodeMatrix << -1,-1,-1,
                    0,-1,-1,
                    1,-1,-1,
                    -1,0,-1,
                    1,0,-1,
                    -1,1,-1,
                    0,1,-1,
                    1,1,-1,
                    -1,-1,0,
                    1,-1,0,
                    -1,1,0,
                    1,1,0,
                    -1,-1,1,
                    0,-1,1,
                    1,-1,1,
                    -1,0,1,
                    1,0,1,
                    -1,1,1,
                    0,1,1,
                    1,1,1;
    }
    else if (dimension == 4)
    {
      if (degree == 1)
      {
        nodeMatrix << -1,-1,-1,-1,
                     1,-1,-1,-1,
                     -1,1,-1,-1,
                     1,1,-1,-1,
                     -1,-1,1,-1,
                     1,-1,1,-1,
                     -1,1,1,-1,
                     1,1,1,-1,
                     -1,-1,-1,1,
                     1,-1,-1,1,
                     -1,1,-1,1,
                     1,1,-1,1,
                     -1,-1,1,1,
                     1,-1,1,1,
                     -1,1,1,1,
                     1,1,1,1;
      }
      else if (degree == 2)
      {
        nodeMatrix << -1,-1,-1,-1,
                    0,-1,-1,-1,
                    1,-1,-1,-1,
                    -1,0,-1,-1,
                    1,0,-1,-1,
                    -1,1,-1,-1,
                    0,1,-1,-1,
                    1,1,-1,-1,
                    -1,-1,0,-1,
                    1,-1,0,-1,
                    -1,1,0,-1,
                    1,1,0,-1,
                    -1,-1,1,-1,
                    0,-1,1,-1,
                    1,-1,1,-1,
                    -1,0,1,-1,
                    1,0,1,-1,
                    -1,1,1,-1,
                    0,1,1,-1,
                    1,1,1,-1,
                    -1,-1,-1,0,
                    1,-1,-1,0,
                    -1,1,-1,0,
                    1,1,-1,0,
                    -1,-1,1,0,
                    1,-1,1,0,
                    -1,1,1,0,
                    1,1,1,0,
                    -1,-1,-1,1,
                    0,-1,-1,1,
                    1,-1,-1,1,
                    -1,0,-1,1,
                    1,0,-1,1,
                    -1,1,-1,1,
                    0,1,-1,1,
                    1,1,-1,1,
                    -1,-1,0,1,
                    1,-1,0,1,
                    -1,1,0,1,
                    1,1,0,1,
                    -1,-1,1,1,
                    0,-1,1,1,
                    1,-1,1,1,
                    -1,0,1,1,
                    1,0,1,1,
                    -1,1,1,1,
                    0,1,1,1,
                    1,1,1,1;
      }
    }
    else if (dimension == 5)
    {
      if (degree == 1)
      {
        nodeMatrix << -1,-1,-1,-1,-1,
                     1,-1,-1,-1,-1,
                     -1,1,-1,-1,-1,
                     1,1,-1,-1,-1,
                     -1,-1,1,-1,-1,
                     1,-1,1,-1,-1,
                     -1,1,1,-1,-1,
                     1,1,1,-1,-1,
                     -1,-1,-1,1,-1,
                     1,-1,-1,1,-1,
                     -1,1,-1,1,-1,
                     1,1,-1,1,-1,
                     -1,-1,1,1,-1,
                     1,-1,1,1,-1,
                     -1,1,1,1,-1,
                     1,1,1,1,-1,
                     -1,-1,-1,-1,1,
                     1,-1,-1,-1,1,
                     -1,1,-1,-1,1,
                     1,1,-1,-1,1,
                     -1,-1,1,-1,1,
                     1,-1,1,-1,1,
                     -1,1,1,-1,1,
                     1,1,1,-1,1,
                     -1,-1,-1,1,1,
                     1,-1,-1,1,1,
                     -1,1,-1,1,1,
                     1,1,-1,1,1,
                     -1,-1,1,1,1,
                     1,-1,1,1,1,
                     -1,1,1,1,1,
                     1,1,1,1,1;
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::setupBasisMatrix(Eigen::MatrixXi& basisMatrix, int degree)
  {
    int dataIndex = 0;
    int shape[NDIM];
    int idx[NDIM];

    for(int dimIndex = 0; dimIndex < NDIM; dimIndex++)
      shape[dimIndex] = degree + 1;

    Lucee::Region<NDIM, int> polyRegion(shape);
    Lucee::RowMajorSequencer<NDIM> polySeq = RowMajorSequencer<NDIM>(polyRegion);

    while(polySeq.step())
    {
      polySeq.fillWithIndex(idx);

      // Compute superlinear degree of this monimial
      int superDegree = 0;
      for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
        superDegree += (idx[dimIndex] > 1)*idx[dimIndex];
      
      if (superDegree < degree + 1)
      {
        for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
          basisMatrix(dataIndex, dimIndex) = idx[dimIndex];
        dataIndex++;
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::computeBasisFunctions(std::vector<blitz::Array<double,NDIM> >& functionVector,
      const Eigen::MatrixXd& nodeList, int degree)
  {
    int numNodes = getSerendipityDimension(degree, NDIM);
    // Matrix to represent basis monomials
    Eigen::MatrixXi basisList = Eigen::MatrixXi::Zero(numNodes, NDIM);
    // Populate basisList according to polyOrder
    setupBasisMatrix(basisList, degree);
    // Compute coefficients for basis functions
    MatrixXd coeffMatrix(numNodes, numNodes);

    // Compute monomial terms evaluated at each nodal coordinate
    for (int nodeIndex = 0; nodeIndex < nodeList.rows(); nodeIndex++)
    {
      for (int basisIndex = 0; basisIndex < basisList.rows(); basisIndex++)
      {
        double coeffProduct = 1.0;
        
        for (int dimIndex = 0; dimIndex < basisList.cols(); dimIndex++)
          coeffProduct = coeffProduct*pow(nodeList(nodeIndex, dimIndex), basisList(basisIndex, dimIndex));
        
        coeffMatrix(nodeIndex,basisIndex) = coeffProduct;
      }
    }

    // Each column of invCoeffMatrix will be coefficient of basis monomial in basisList
    MatrixXd invCoeffMatrix = coeffMatrix.inverse();

    blitz::TinyVector<int, NDIM> polynomialShape(maxPower);
    blitz::TinyVector<int, NDIM> polynomialCoord;
    blitz::Array<double, NDIM> polynomialArray(polynomialShape);

    // Put together representations of the basis functions into functionVector
    // Loop over the coefficient matrix (each column represents a basis function)
    for (int basisIndex = 0; basisIndex < invCoeffMatrix.cols(); basisIndex++)
    {
      // Reset the polynomialArray
      polynomialArray = 0;
      // Loop over the monomial terms that make up the basis function
      for (int monomialIndex = 0; monomialIndex < invCoeffMatrix.rows(); monomialIndex++)
      {
        polynomialCoord = 0;
        // Store the basis function in polynomialArray by storing its coefficient at
        // a location specified by the monimial term powers
        for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
          polynomialCoord(dimIndex) = basisList(monomialIndex, dimIndex);

        polynomialArray(polynomialCoord) = polynomialArray(polynomialCoord) + invCoeffMatrix(monomialIndex, basisIndex);
      }
      // Store the newly represented basis function into functionVector
      functionVector.push_back(polynomialArray.copy());
    }
  }

  template <unsigned NDIM>
  double
  SerendipityElement<NDIM>::evalPolynomial(const blitz::Array<double,NDIM>& polyCoeffs, const VectorXd& nodeCoords) const
  {
    double totalSum = 0.0;
    double monomialTerm;
    int shape[NDIM];
    int idx[NDIM];
    blitz::TinyVector<int, NDIM> polyCoord;

    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
      shape[dimIndex] = polyCoeffs.extent(dimIndex);

    Lucee::Region<NDIM, int> polyRegion(shape);
    Lucee::RowMajorSequencer<NDIM> polySeq = RowMajorSequencer<NDIM>(polyRegion);
    // Loop over each element of polyCoeffs
    while(polySeq.step())
    {
      polySeq.fillWithIndex(idx);
      for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
        polyCoord(dimIndex) = idx[dimIndex];
      
      if (polyCoeffs(polyCoord) != 0.0)
      {
        monomialTerm = polyCoeffs(polyCoord);
        // Compute monomial term associated with polyCoord
        for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
          monomialTerm *= pow(nodeCoords(dimIndex), idx[dimIndex]);
        totalSum += monomialTerm;
      }
    }

    return totalSum;
  }

  template <unsigned NDIM>
  blitz::Array<double, NDIM>
  SerendipityElement<NDIM>::computePolynomialDerivative(const blitz::Array<double,NDIM>& poly, int dir)
  {
    int shape[NDIM];
    int idx[NDIM];
    blitz::TinyVector<int, NDIM> polyCoord;
    blitz::Array<double, NDIM> polyResult(poly.shape());
    polyResult = 0;

    for(int dimIndex = 0; dimIndex < NDIM; dimIndex++)
      shape[dimIndex] = poly.extent(dimIndex);

    Lucee::Region<NDIM, int> polyRegion(shape);

    Lucee::RowMajorSequencer<NDIM> polySeq = RowMajorSequencer<NDIM>(polyRegion);
    // Loop over each element of polyCoeffs
    while(polySeq.step())
    {
      polySeq.fillWithIndex(idx);
      for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
        polyCoord(dimIndex) = idx[dimIndex];
      
      if (polyCoord(dir) != 0 && poly(polyCoord) != 0.0)
      {
        double resultCoeff = polyCoord(dir)*poly(polyCoord);
        // Compute coordinate of derivative term
        polyCoord(dir) = polyCoord(dir) - 1;
        // Store new coefficient
        polyResult(polyCoord) += resultCoeff;
      }
    }

    return polyResult;
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::computeMass(Eigen::MatrixXd& resultMatrix)
  {
    for (int kIndex = 0; kIndex < resultMatrix.rows(); kIndex++)
    {
      for (int mIndex = 0; mIndex < resultMatrix.cols(); mIndex++)
      {
        // Reset integration result
        double integrationResult = 0.0;

        // Loop over 3d gauss points to evaluate integral using gaussian quadrature
        // result = weight(node) * f1(node) * f2(node)
        for (int gaussIndex = 0; gaussIndex < gaussNodeList.rows(); gaussIndex++)
          integrationResult += gaussNodeList(gaussIndex, NDIM)*functionEvaluations(gaussIndex, kIndex)*functionEvaluations(gaussIndex, mIndex);
        
        resultMatrix(kIndex, mIndex) = integrationResult;
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::computeFaceMass(int dir, Eigen::MatrixXd& lowerResultMatrix, Eigen::MatrixXd& upperResultMatrix)
  {
    std::vector<int> surfLowerNodeNums(lowerResultMatrix.cols(),0);
    std::vector<int> surfUpperNodeNums(upperResultMatrix.cols(),0);

    getSurfLowerNodeNums(dir,surfLowerNodeNums);
    getSurfUpperNodeNums(dir,surfUpperNodeNums);

    for (int kIndex = 0; kIndex < upperResultMatrix.rows(); kIndex++)
    {
      // Note that mIndex will be those of the ones in lower/upper node nums
      for (int mIndex = 0; mIndex < upperResultMatrix.cols(); mIndex++)
      {
        int lowerNodeNum = surfLowerNodeNums[mIndex];
        int upperNodeNum = surfUpperNodeNums[mIndex];
        double integrationResultL = 0.0;
        double integrationResultU = 0.0;
        
        for (int testIndex = 0; testIndex < gaussNodeListUpperSurf[dir].rows(); testIndex++)
        {
          integrationResultU += upperSurfaceEvaluations[dir](testIndex, upperNodeNum)*
            upperSurfaceEvaluations[dir](testIndex, kIndex)*gaussNodeListUpperSurf[dir](testIndex, NDIM);
          integrationResultL += lowerSurfaceEvaluations[dir](testIndex,lowerNodeNum)*
            lowerSurfaceEvaluations[dir](testIndex, kIndex)*gaussNodeListLowerSurf[dir](testIndex, NDIM);
        }

        upperResultMatrix(kIndex,mIndex) = integrationResultU;
        lowerResultMatrix(kIndex,mIndex) = integrationResultL;
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::computeStiffness(const blitz::Array<double, 3>& functionDerivative,
    Eigen::MatrixXd& resultMatrix)
  {
    for (int kIndex = 0; kIndex < resultMatrix.rows(); kIndex++)
    {
      for (int mIndex = 0; mIndex < resultMatrix.cols(); mIndex++)
      {
        // Reset integration result
        double integrationResult = 0.0;
        
        // Loop over 3d gauss points to evaluate integral using gaussian quadrature
        for (int gaussIndex = 0; gaussIndex < gaussNodeList.rows(); gaussIndex++)
          for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
            integrationResult += gaussNodeList(gaussIndex,NDIM)*functionDerivative(gaussIndex, kIndex, dimIndex)*functionDerivative(gaussIndex, mIndex, dimIndex)*4.0/dq2[dimIndex];
        
        resultMatrix(kIndex,mIndex) = integrationResult;
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::computeGradStiffness(const blitz::Array<double, 3>& functionDerivative,
    int dir, Eigen::MatrixXd& resultMatrix)
  {
    for (int kIndex = 0; kIndex < resultMatrix.rows(); kIndex++)
    {
      for (int mIndex = 0; mIndex < resultMatrix.cols(); mIndex++)
      {
        // Reset integration result
        double integrationResult = 0.0;

        // Loop over 3d gauss points to evaluate integral using gaussian quadrature
        for (int gaussIndex = 0; gaussIndex < gaussNodeList.rows(); gaussIndex++)
          integrationResult += gaussNodeList(gaussIndex,NDIM)*functionDerivative(gaussIndex, kIndex, dir)*functionEvaluations(gaussIndex, mIndex);

        resultMatrix(kIndex, mIndex) = integrationResult;
      }
    }
  }
 
  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::setupMomentMatrices()
  {
    unsigned nodesPerSide = getNumSurfLowerNodes(1);
    std::vector<int> surfLowerNodeNums(nodesPerSide,1);
    getSurfLowerNodeNums(1,surfLowerNodeNums);

    maxMoment = 3;
    momMatrix = std::vector<Eigen::MatrixXd>(maxMoment+1);

    // Initialize Matrices
    for (int momentVal = 0; momentVal < maxMoment + 1; momentVal++)
      momMatrix[momentVal] = Eigen::MatrixXd::Zero(nodesPerSide, functionEvaluations.cols());

    // Evaluate Integral y^p*phi(x)psi(x,y)dA
    for (int j = 0; j < nodesPerSide; j++)
    {
      // 1-d basis function in this loop
      int lowerNodeNum = surfLowerNodeNums[j];
      for (int k = 0; k < functionEvaluations.cols(); k++)
      {
        for (int nodeIndex = 0; nodeIndex < gaussNodeList.rows(); nodeIndex++)
        {
          double gaussWeight  = gaussNodeList(nodeIndex,NDIM);
          double yCoord       = gaussNodeList(nodeIndex,1);
          unsigned rollingIndex = nodeIndex % numGaussPoints;

          double yPowerFactor = 1.0;
          for (int momentVal = 0; momentVal <= maxMoment; momentVal++)
          {
            if (momentVal != 0)
              yPowerFactor *= yCoord;
            momMatrix[momentVal](j,k) += gaussWeight*yPowerFactor*functionEvaluations(nodeIndex, k)*
              lowerSurfaceEvaluations[1](rollingIndex, lowerNodeNum);
          }
        }
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::computeTransformationScales()
  {
    // Temporary code for Jacobian calculation:
    Eigen::MatrixXd refVertexList = Eigen::MatrixXd::Zero(getSerendipityDimension(1,NDIM),NDIM);
    // Get vertices (nodes of polyOrder 1 element)
    getNodeList(refVertexList, 1, NDIM);
    std::vector<blitz::Array<double,NDIM> > vertexFunctionVector;
    computeBasisFunctions(vertexFunctionVector, refVertexList, 1);

    Eigen::MatrixXd physicalVertexList = refVertexList;
    
    std::vector<blitz::Array<double,NDIM> > transformationVector;

    blitz::TinyVector<int, NDIM> polynomialShape(maxPower);
    blitz::Array<double, NDIM> polynomialArray(polynomialShape);

    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      // Reset polynomialArray
      polynomialArray = 0;
      for (int nodeIndex = 0; nodeIndex < physicalVertexList.rows(); nodeIndex++)
      {
        polynomialArray = polynomialArray + physicalVertexList(nodeIndex, dimIndex)*vertexFunctionVector[nodeIndex];
      }

      transformationVector.push_back(polynomialArray.copy());
    }

    std::vector<blitz::Array<double,NDIM> > jacobianVector;
    Eigen::MatrixXd volQuadJacobian(NDIM*NDIM, gaussNodeList.rows());
    //Eigen::MatrixXd surfQuadJacobian(NDIM*NDIM, gaussNodeList.rows());
    std::vector<Eigen::MatrixXd> surfUpperJacobianVector(NDIM);
    std::vector<Eigen::MatrixXd> surfLowerJacobianVector(NDIM);
    
    // Allocate Eigen matrices
    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      surfUpperJacobianVector[dimIndex] = Eigen::MatrixXd(NDIM*NDIM, gaussNodeListUpperSurf[dimIndex].rows());
      surfLowerJacobianVector[dimIndex] = Eigen::MatrixXd(NDIM*NDIM, gaussNodeListLowerSurf[dimIndex].rows());
    }

    // Evaluate derivatives of basis functions at various vol and surf quadrature points
    for (int rowIndex = 0; rowIndex < NDIM; rowIndex++)
    {
      for (int colIndex = 0; colIndex < NDIM; colIndex++)
      {
        // Matrix Element (rowIndex,colIndex) is rowIndex reference coordinate derivative on
        // colIndex physical coordinate
        blitz::Array<double, NDIM> polyDeriv = computePolynomialDerivative(transformationVector[colIndex], rowIndex);
        for (int volIndex = 0; volIndex < gaussNodeList.rows(); volIndex++)
        {
          Eigen::VectorXd nodeVec(NDIM);
          for(int i = 0; i < NDIM; i++)
            nodeVec(i) = gaussNodeList(volIndex, i);
          // Evaluate and store basis function derivative at quadrature point
          volQuadJacobian(rowIndex*NDIM + colIndex, volIndex) = evalPolynomial(polyDeriv, nodeVec);
        }

        // Jacobian matrices evaluated at surface quadrature points
        for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
        {
          for (int surfIndex = 0; surfIndex < gaussNodeListUpperSurf[dimIndex].rows(); surfIndex++)
          {
            Eigen::VectorXd nodeVec(NDIM);
            for(int i = 0; i < NDIM; i++)
              nodeVec(i) = gaussNodeListUpperSurf[dimIndex](surfIndex, i);

            surfUpperJacobianVector[dimIndex](rowIndex*NDIM + colIndex, surfIndex) = evalPolynomial(polyDeriv, nodeVec);

            for(int i = 0; i < NDIM; i++)
              nodeVec(i) = gaussNodeListLowerSurf[dimIndex](surfIndex, i);

            surfLowerJacobianVector[dimIndex](rowIndex*NDIM + colIndex, surfIndex) = evalPolynomial(polyDeriv, nodeVec);
          }
        }
      }
    }

    // Compute jacobians of various matrices for volume integrals
    for (int volIndex = 0; volIndex < gaussNodeList.rows(); volIndex++)
    {
      Eigen::MatrixXd jacobianMatrix(NDIM, NDIM);
      for (int i = 0; i < NDIM; i++)
        for (int j = 0; j < NDIM; j++)
          jacobianMatrix(i, j) = volQuadJacobian(i*NDIM + j, volIndex);

      // Multiply corresponding quadrature weight by det(J)
      gaussNodeList(volIndex, NDIM) *= jacobianMatrix.determinant();

      // Store jacobian or inverse of jacobian matrix for use in computing stiffness matrices
    }

    // Compute jacobians of various matrices for surface integrals
    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      for (int surfIndex = 0; surfIndex < gaussNodeListUpperSurf[dimIndex].rows(); surfIndex++)
      {
        Eigen::MatrixXd jacobianMatrix(NDIM, NDIM);
        for (int i = 0; i < NDIM; i++)
          for (int j = 0; j < NDIM; j++)
            jacobianMatrix(i, j) = surfUpperJacobianVector[dimIndex](i*NDIM + j, surfIndex);

        // Multiply corresponding quadrature weight by det(J)
        gaussNodeListUpperSurf[dimIndex](surfIndex, NDIM) *= jacobianMatrix.determinant();

        // Store jacobian or inverse of jacobian matrix for use in computing stiffness matrices
        
        // Do the same for the lower surface
        for (int i = 0; i < NDIM; i++)
          for (int j = 0; j < NDIM; j++)
            jacobianMatrix(i, j) = surfLowerJacobianVector[dimIndex](i*NDIM + j, surfIndex);
        gaussNodeListLowerSurf[dimIndex](surfIndex, NDIM) *= jacobianMatrix.determinant();
      }
    }
  }

  template <unsigned NDIM>
  int
  SerendipityElement<NDIM>::getSerendipityDimension(int degree, int dimension) const
  {
    int upperBound = std::min(dimension, (int) std::floor(degree/2.0));
    int totalTerms = 0;

    for (int d = 0; d < upperBound + 1; d++)
    {
      totalTerms = totalTerms + pow(2,dimension-d)*factorial(dimension)/(factorial(dimension-d)*factorial(d))*
        factorial(degree-d)/(factorial(degree-2*d)*factorial(d));
    }

    return totalTerms;
  }

  template <unsigned NDIM>
  int
  SerendipityElement<NDIM>::factorial(int n) const
  {
    return (n == 1 || n == 0) ? 1 : factorial(n-1)*n;
  }
  
  // instantiations
  template class SerendipityElement<1>;
  template class SerendipityElement<2>;
  template class SerendipityElement<3>;
  template class SerendipityElement<4>;
  template class SerendipityElement<5>;
}
