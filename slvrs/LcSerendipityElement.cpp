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
    maxPower = 3*polyOrder;
    
    if (NDIM == 2)
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
    else if (NDIM == 3)
    {
      if (polyOrder < 5)
        this->setNumNodes(getSerendipityDimension(polyOrder, NDIM));
      else
      {
        Lucee::Except lce("SerendipityElement: Degree must be 1, 2, 3, or 4.");
        lce << " Provided " << polyOrder << " instead";
        throw lce;
      }
      setupMatrices();
    }
    else
    {
      Lucee::Except lce("SerendipityElement: NDIM must be 2 or 3.");
      lce << " Provided " << NDIM << " instead";
      throw lce;
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getExclusiveNodeIndices(std::vector<int>& ndIds)
  {
    if (NDIM == 2)
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
      Lucee::Except lce("SerendipityElement::getExclusiveNodeIndices NDIM must be 2.");
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
    throw Lucee::Except("SerendipityElement::getNumGlobalNodes: Not implemented!");
    return 0;
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getLocalToGlobal(std::vector<int>& lgMap) const
  {
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
    unsigned nodesPerRow2D[nodesPerSide];
    unsigned globalNodesPerRow2D[nodesPerSide];
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

    unsigned nodeIndices[nodesPerCell2D];

    throw Lucee::Except("SerendipityElement::getLocalToGlobal: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getSurfLowerLocalToGlobal(unsigned dim,
    std::vector<int>& lgMap) const
  {
    throw Lucee::Except("SerendipityElement::getSurfLowerLocalToGlobal: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getSurfUpperLocalToGlobal(unsigned dim,
    std::vector<int>& lgMap) const
  {
    throw Lucee::Except("SerendipityElement::getSurfUpperLocalToGlobal: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getSurfLowerNodeNums(unsigned dir,
    std::vector<int>& nodeNum) const
  {
    if (NDIM == 2)
    {
      if (polyOrder == 1)
      {
        if (dir == 0)
        {
          nodeNum[0] = 0;
          nodeNum[1] = 3;
        }
        else if (dir == 1)
        {
          nodeNum[0] = 0;
          nodeNum[1] = 1;
        }
      }
      else if (polyOrder == 2)
      {
        if (dir == 0)
        {
          nodeNum[0] = 0;
          nodeNum[1] = 7;
          nodeNum[2] = 6;
        }
        else if (dir == 1)
        {
          nodeNum[0] = 0;
          nodeNum[1] = 1;
          nodeNum[2] = 2;
        }
      }
      else if (polyOrder == 3)
      {
        if (dir == 0)
        {
          nodeNum[0] = 0;
          nodeNum[1] = 11;
          nodeNum[2] = 10;
          nodeNum[3] = 9;
        }
        else if (dir == 1)
        {
          nodeNum[0] = 0;
          nodeNum[1] = 1;
          nodeNum[2] = 2;
          nodeNum[3] = 3;
        }
      }
    }
    else if (NDIM == 3)
    {
      if (polyOrder == 1)
      {
        if (dir == 0)
        {
          nodeNum[0] = 0;
          nodeNum[1] = 3;
          nodeNum[2] = 4;
          nodeNum[3] = 7;
        }
        else if (dir == 1)
        {
          nodeNum[0] = 0;
          nodeNum[1] = 1;
          nodeNum[2] = 4;
          nodeNum[3] = 5;
        }
        else if (dir == 2)
        {
          nodeNum[0] = 0;
          nodeNum[1] = 1;
          nodeNum[2] = 3;
          nodeNum[3] = 2;
        }
      }
      else if (polyOrder == 2)
      {
        if (dir == 0)
        {
          nodeNum[0] = 0;
          nodeNum[1] = 7;
          nodeNum[2] = 6;
          nodeNum[3] = 8;
          nodeNum[4] = 11;
          nodeNum[5] = 12;
          nodeNum[6] = 19;
          nodeNum[7] = 18;
        }
        else if (dir == 1)
        {
          nodeNum[0] = 0;
          nodeNum[1] = 1;
          nodeNum[2] = 2;
          nodeNum[3] = 8;
          nodeNum[4] = 9;
          nodeNum[5] = 12;
          nodeNum[6] = 13;
          nodeNum[7] = 14;
        }
        else if (dir == 2)
        {
          nodeNum[0] = 0;
          nodeNum[1] = 1;
          nodeNum[2] = 2;
          nodeNum[3] = 7;
          nodeNum[4] = 3;
          nodeNum[5] = 6;
          nodeNum[6] = 5;
          nodeNum[7] = 4;
        }
      }
      else if (polyOrder == 3)
      {
        if (dir == 0)
        {
          nodeNum[0] = 0;
          nodeNum[1] = 11;
          nodeNum[2] = 10;
          nodeNum[3] = 9;
          nodeNum[4] = 12;
          nodeNum[5] = 15;
          nodeNum[6] = 16;
          nodeNum[7] = 19;
          nodeNum[8] = 20;
          nodeNum[9] = 31;
          nodeNum[10] = 30;
          nodeNum[11] = 29;
        }
        else if (dir == 1)
        {
          nodeNum[0] = 0;
          nodeNum[1] = 1;
          nodeNum[2] = 2;
          nodeNum[3] = 3;
          nodeNum[4] = 12;
          nodeNum[5] = 13;
          nodeNum[6] = 16;
          nodeNum[7] = 17;
          nodeNum[8] = 20;
          nodeNum[9] = 21;
          nodeNum[10] = 22;
          nodeNum[11] = 23;
        }
        else if (dir == 2)
        {
          nodeNum[0] = 0;
          nodeNum[1] = 1;
          nodeNum[2] = 2;
          nodeNum[3] = 3;
          nodeNum[4] = 11;
          nodeNum[5] = 4;
          nodeNum[6] = 10;
          nodeNum[7] = 5;
          nodeNum[8] = 9;
          nodeNum[9] = 8;
          nodeNum[10] = 7;
          nodeNum[11] = 6;
        }
      }
      else if (polyOrder == 4)
      {
        if (dir == 0)
        {
          nodeNum[0] = 0;
          nodeNum[1] = 15;
          nodeNum[2] = 14;
          nodeNum[3] = 13;
          nodeNum[4] = 12;
          nodeNum[5] = 17;
          nodeNum[6] = 20;
          nodeNum[7] = 21;
          nodeNum[8] = 28;
          nodeNum[9] = 27;
          nodeNum[10] = 29;
          nodeNum[11] = 32;
          nodeNum[12] = 33;
          nodeNum[13] = 48;
          nodeNum[14] = 47;
          nodeNum[15] = 46;
          nodeNum[16] = 45;
        }
        else if (dir == 1)
        {
          nodeNum[0] = 0;
          nodeNum[1] = 1;
          nodeNum[2] = 2;
          nodeNum[3] = 3;
          nodeNum[4] = 4;
          nodeNum[5] = 17;
          nodeNum[6] = 18;
          nodeNum[7] = 21;
          nodeNum[8] = 22;
          nodeNum[9] = 23;
          nodeNum[10] = 29;
          nodeNum[11] = 30;
          nodeNum[12] = 33;
          nodeNum[13] = 34;
          nodeNum[14] = 35;
          nodeNum[15] = 36;
          nodeNum[16] = 37;
        }
        else if (dir == 2)
        {
          nodeNum[0] = 0;
          nodeNum[1] = 1;
          nodeNum[2] = 2;
          nodeNum[3] = 3;
          nodeNum[4] = 4;
          nodeNum[5] = 15;
          nodeNum[6] = 5;
          nodeNum[7] = 14;
          nodeNum[8] = 16;
          nodeNum[9] = 6;
          nodeNum[10] = 13;
          nodeNum[11] = 7;
          nodeNum[12] = 12;
          nodeNum[13] = 11;
          nodeNum[14] = 10;
          nodeNum[15] = 9;
          nodeNum[16] = 8;
        }
      }

    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getSurfUpperNodeNums(unsigned dir,
    std::vector<int>& nodeNum) const
  {
    if (NDIM == 2)
    {
      if (polyOrder == 1)
      {
        if (dir == 0)
        {
          nodeNum[0] = 1;
          nodeNum[1] = 2;
        }
        else if (dir == 1)
        {
          nodeNum[0] = 3;
          nodeNum[1] = 2;
        }
      }
      else if (polyOrder == 2)
      {
        if (dir == 0)
        {
          nodeNum[0] = 2;
          nodeNum[1] = 3;
          nodeNum[2] = 4;
        }
        else if (dir == 1)
        {
          nodeNum[0] = 6;
          nodeNum[1] = 5;
          nodeNum[2] = 4;
        }
      }
      else if (polyOrder == 3)
      {
        if (dir == 0)
        {
          nodeNum[0] = 3;
          nodeNum[1] = 4;
          nodeNum[2] = 5;
          nodeNum[3] = 6;
        }
        else if (dir == 1)
        {
          nodeNum[0] = 9;
          nodeNum[1] = 8;
          nodeNum[2] = 7;
          nodeNum[3] = 6;
        }
      }
    }
    else if (NDIM == 3)
    {
      if (polyOrder == 1)
      {
        if (dir == 0)
        {
          nodeNum[0] = 1;
          nodeNum[1] = 2;
          nodeNum[2] = 5;
          nodeNum[3] = 6;
        }
        else if (dir == 1)
        {
          nodeNum[0] = 3;
          nodeNum[1] = 2;
          nodeNum[2] = 7;
          nodeNum[3] = 6;
        }
        else if (dir == 2)
        {
          nodeNum[0] = 4;
          nodeNum[1] = 5;
          nodeNum[2] = 7;
          nodeNum[3] = 6;
        }
      }
      else if (polyOrder == 2)
      {
        if (dir == 0)
        {
          nodeNum[0] = 2;
          nodeNum[1] = 3;
          nodeNum[2] = 4;
          nodeNum[3] = 9;
          nodeNum[4] = 10;
          nodeNum[5] = 14;
          nodeNum[6] = 15;
          nodeNum[7] = 16;
        }
        else if (dir == 1)
        {
          nodeNum[0] = 6;
          nodeNum[1] = 5;
          nodeNum[2] = 4;
          nodeNum[3] = 11;
          nodeNum[4] = 10;
          nodeNum[5] = 18;
          nodeNum[6] = 17;
          nodeNum[7] = 16;
        }
        else if (dir == 2)
        {
          nodeNum[0] = 12;
          nodeNum[1] = 13;
          nodeNum[2] = 14;
          nodeNum[3] = 19;
          nodeNum[4] = 15;
          nodeNum[5] = 18;
          nodeNum[6] = 17;
          nodeNum[7] = 16;
        }
      }
      else if (polyOrder == 3)
      {
        if (dir == 0)
        {
          nodeNum[0] = 3;
          nodeNum[1] = 4;
          nodeNum[2] = 5;
          nodeNum[3] = 6;
          nodeNum[4] = 13;
          nodeNum[5] = 14;
          nodeNum[6] = 17;
          nodeNum[7] = 18;
          nodeNum[8] = 23;
          nodeNum[9] = 24;
          nodeNum[10] = 25;
          nodeNum[11] = 26;
        }
        else if (dir == 1)
        {
          nodeNum[0] = 9;
          nodeNum[1] = 8;
          nodeNum[2] = 7;
          nodeNum[3] = 6;
          nodeNum[4] = 15;
          nodeNum[5] = 14;
          nodeNum[6] = 19;
          nodeNum[7] = 18;
          nodeNum[8] = 29;
          nodeNum[9] = 28;
          nodeNum[10] = 27;
          nodeNum[11] = 26;
        }
        else if (dir == 2)
        {
          nodeNum[0] = 20;
          nodeNum[1] = 21;
          nodeNum[2] = 22;
          nodeNum[3] = 23;
          nodeNum[4] = 31;
          nodeNum[5] = 24;
          nodeNum[6] = 30;
          nodeNum[7] = 25;
          nodeNum[8] = 29;
          nodeNum[9] = 28;
          nodeNum[10] = 27;
          nodeNum[11] = 26;
        }
      }
      else if (polyOrder == 4)
      {
        if (dir == 0)
        {
          nodeNum[0] = 4;
          nodeNum[1] = 5;
          nodeNum[2] = 6;
          nodeNum[3] = 7;
          nodeNum[4] = 8;
          nodeNum[5] = 18;
          nodeNum[6] = 19;
          nodeNum[7] = 23;
          nodeNum[8] = 24;
          nodeNum[9] = 25;
          nodeNum[10] = 30;
          nodeNum[11] = 31;
          nodeNum[12] = 37;
          nodeNum[13] = 38;
          nodeNum[14] = 39;
          nodeNum[15] = 40;
          nodeNum[16] = 41;
        }
        else if (dir == 1)
        {
          nodeNum[0] = 12;
          nodeNum[1] = 11;
          nodeNum[2] = 10;
          nodeNum[3] = 9;
          nodeNum[4] = 8;
          nodeNum[5] = 20;
          nodeNum[6] = 19;
          nodeNum[7] = 27;
          nodeNum[8] = 26;
          nodeNum[9] = 25;
          nodeNum[10] = 32;
          nodeNum[11] = 31;
          nodeNum[12] = 45;
          nodeNum[13] = 44;
          nodeNum[14] = 43;
          nodeNum[15] = 42;
          nodeNum[16] = 41;
        }
        else if (dir == 2)
        {
          nodeNum[0] = 33;
          nodeNum[1] = 34;
          nodeNum[2] = 35;
          nodeNum[3] = 36;
          nodeNum[4] = 37;
          nodeNum[5] = 48;
          nodeNum[6] = 38;
          nodeNum[7] = 47;
          nodeNum[8] = 49;
          nodeNum[9] = 39;
          nodeNum[10] = 46;
          nodeNum[11] = 40;
          nodeNum[12] = 45;
          nodeNum[13] = 44;
          nodeNum[14] = 43;
          nodeNum[15] = 42;
          nodeNum[16] = 41;
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
    double xc[3];

    grid.getCentroid(xc);

    // Loop over all node locations on reference element and convert them
    // to appropriate coordinates
    for (int i = 0; i < this->getNumNodes(); i++)
      for (int dim = 0; dim < NDIM; dim++)
        nodeCoords(i, dim) = xc[dim] + nodeList(i, dim)*0.5*dq[dim];
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
      if (NDIM == 2)
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
      if (NDIM == 2)
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
      if (NDIM == 2)
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
    if (p > maxMoment)
    {
      // moments higher than 3 not supported for now
      Lucee::Except lce("SerendipityElement::getMomentMatrix: Moment matrix of order ");
      lce << p << " not supported";
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
    getNodeList(nodeList);

    std::vector<blitz::Array<double,NDIM> > functionVector;
    computeBasisFunctions(functionVector);

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
    Eigen::MatrixXd gaussNodeListSurf = Eigen::MatrixXd::Zero(totalSurfaceGaussNodes, NDIM);
    
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
      #include <LcSerendipityElementFaceToInteriorOutput>
    }

    computeMass(refMass);
    computeStiffness(functionDEvaluations,refStiffness);

    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      computeFaceMass(functionVector, dimIndex, refFaceMassLower[dimIndex], refFaceMassUpper[dimIndex]);
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
  SerendipityElement<NDIM>::getNodeList(Eigen::MatrixXd& nodeMatrix)
  {
    if (NDIM == 2)
    {
      if (polyOrder == 1)
      {
        nodeMatrix << -1,-1,
                      1,-1,
                      1,1,
                      -1,1;
      }
      else if (polyOrder == 2)
      {
        nodeMatrix << -1,-1,
                    0,-1,
                    1,-1,
                    1,0,
                    1,1,
                    0,1,
                    -1,1,
                    -1,0;
      }
      else if (polyOrder == 3)
      {
        nodeMatrix << -1,-1,
                    -1/3.0,-1,
                    1/3.0,-1,
                    1,-1,
                    1,-1/3.0,
                    1,1/3.0,
                    1,1,
                    1/3.0,1,
                    -1/3.0,1,
                    -1,1,
                    -1,1/3.0,
                    -1,-1/3.0;
      }
    }
    else if (NDIM == 3)
    {
      if (polyOrder == 1) {
        nodeMatrix << -1,-1,-1,
                     1,-1,-1,
                     1,1,-1,
                     -1,1,-1,
                     -1,-1,1,
                     1,-1,1,
                     1,1,1,
                     -1,1,1;
      }
      else if (polyOrder == 2)
      {
        nodeMatrix << -1,-1,-1,
                    0,-1,-1,
                    1,-1,-1,
                    1,0,-1,
                    1,1,-1,
                    0,1,-1,
                    -1,1,-1,
                    -1,0,-1,// End of Bottom Layer
                    -1,-1,0,
                    1,-1,0,
                    1,1,0,
                    -1,1,0, // End of Middle Layer
                    -1,-1,1,
                    0,-1,1,
                    1,-1,1,
                    1,0,1,
                    1,1,1,
                    0,1,1,
                    -1,1,1,
                    -1,0,1; // End of Top Layer
      }
      else if (polyOrder == 3)
      {
        nodeMatrix << -1,-1,-1,
                    -1/3.0,-1,-1,
                    1/3.0,-1,-1,
                    1,-1,-1,
                    1,-1/3.0,-1,
                    1,1/3.0,-1,
                    1,1,-1,
                    1/3.0,1,-1,
                    -1/3.0,1,-1,
                    -1,1,-1,
                    -1,1/3.0,-1,
                    -1,-1/3.0,-1,// End of Bottom Layer
                    -1,-1,-1/3.0,
                    1,-1,-1/3.0,
                    1,1,-1/3.0,
                    -1,1,-1/3.0, // End of 2nd Layer
                    -1,-1,1/3.0,
                    1,-1,1/3.0,
                    1,1,1/3.0,
                    -1,1,1/3.0,  // End of 3rd Layer
                    -1,-1,1,
                    -1/3.0,-1,1,
                    1/3.0,-1,1,
                    1,-1,1,
                    1,-1/3.0,1,
                    1,1/3.0,1,
                    1,1,1,
                    1/3.0,1,1,
                    -1/3.0,1,1,
                    -1,1,1,
                    -1,1/3.0,1,
                    -1,-1/3.0,1; // End of Top Layer
      }
      else if (polyOrder == 4)
      {
        nodeMatrix << -1,-1,-1,
                      -0.5,-1,-1,
                      0,-1,-1,
                      0.5,-1,-1,
                      1,-1,-1,
                      1,-0.5,-1,
                      1,0,-1,
                      1,0.5,-1,
                      1,1,-1,
                      0.5,1,-1,
                      0,1,-1,
                      -0.5,1,-1,
                      -1,1,-1,
                      -1,0.5,-1,
                      -1,0,-1,
                      -1,-0.5,-1,
                      0,0,-1,
                      -1,-1,-0.5,
                      1,-1,-0.5,
                      1,1,-0.5,
                      -1,1,-0.5,
                      -1,-1,0,
                      0,-1,0,
                      1,-1,0,
                      1,0,0,
                      1,1,0,
                      0,1,0,
                      -1,1,0,
                      -1,0,0,
                      -1,-1,0.5,
                      1,-1,0.5,
                      1,1,0.5,
                      -1,1,0.5,
                      -1,-1,1,
                      -0.5,-1,1,
                      0,-1,1,
                      0.5,-1,1,
                      1,-1,1,
                      1,-0.5,1,
                      1,0,1,
                      1,0.5,1,
                      1,1,1,
                      0.5,1,1,
                      0,1,1,
                      -0.5,1,1,
                      -1,1,1,
                      -1,0.5,1,
                      -1,0,1,
                      -1,-0.5,1,
                      0,0,1;
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::setupBasisMatrix(Eigen::MatrixXi& basisMatrix)
  {
    int dataIndex = 0;
    int shape[NDIM];
    int idx[NDIM];

    for(int dimIndex = 0; dimIndex < NDIM; dimIndex++)
      shape[dimIndex] = polyOrder + 1;

    Lucee::Region<NDIM, int> polyRegion(shape);
    Lucee::RowMajorSequencer<NDIM> polySeq = RowMajorSequencer<NDIM>(polyRegion);

    while(polySeq.step())
    {
      polySeq.fillWithIndex(idx);

      // Compute superlinear degree of this monimial
      int superDegree = 0;
      for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
        superDegree += (idx[dimIndex] > 1)*idx[dimIndex];
      
      if (superDegree < polyOrder + 1)
      {
        for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
          basisMatrix(dataIndex, dimIndex) = idx[dimIndex];
        dataIndex++;
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::computeBasisFunctions(std::vector<blitz::Array<double,NDIM> >& functionVector)
  {
    // Matrix to represent basis monomials
    Eigen::MatrixXi basisList = Eigen::MatrixXi::Zero(this->getNumNodes(), NDIM);
    // Populate basisList according to polyOrder
    setupBasisMatrix(basisList);
    // Compute coefficients for basis functions
    MatrixXd coeffMatrix(this->getNumNodes(), this->getNumNodes());

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
      // Reset the polynomial3DArray
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
  SerendipityElement<NDIM>::evalPolynomial(const blitz::Array<double,NDIM>& polyCoeffs, const VectorXd& nodeCoords)
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
  SerendipityElement<NDIM>::computeFaceMass(const std::vector<blitz::Array<double,NDIM> >& functionVector, int dir,
    Eigen::MatrixXd& lowerResultMatrix, Eigen::MatrixXd& upperResultMatrix)
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
