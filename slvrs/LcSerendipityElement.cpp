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
    // polyOrder is a misnomer for serendipity elements when we really want the degree
    basisDegree = polyOrder;
    
    if (NDIM == 2)
    {
      if (basisDegree == 1)
        this->setNumNodes(4);
      else if (basisDegree == 2)
        this->setNumNodes(8);
      else if (basisDegree == 3)
        this->setNumNodes(12);
      else
      {
        Lucee::Except lce("SerendipityElement: Degree must be 1, 2, or 3.");
        lce << " Provided " << basisDegree << " instead";
        throw lce;
      }
      setupMatrices();
    }
    else if (NDIM == 3)
    {
      if (basisDegree == 1)
        this->setNumNodes(8);
      else if (basisDegree == 2)
        this->setNumNodes(20);
      else if (basisDegree == 3)
        this->setNumNodes(32);
      else if (basisDegree == 4)
        this->setNumNodes(50);
      else
      {
        Lucee::Except lce("SerendipityElement: Degree must be 1, 2, or 3.");
        lce << " Provided " << basisDegree << " instead";
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
    if (NDIM == 2)
    {
      if (basisDegree == 1)
        return 2;
      else if (basisDegree == 2)
        return 3;
      else if (basisDegree == 3)
        return 4;
    }
    else if (NDIM ==3)
    {
      if (basisDegree == 1)
        return 4;
      else if (basisDegree == 2)
        return 8;
      else if (basisDegree == 3)
        return 12;
      else if (basisDegree == 4)
        return 17;
    }
  }

  template <unsigned NDIM>
  unsigned
  SerendipityElement<NDIM>::getNumSurfUpperNodes(unsigned dir) const
  {
    if (NDIM == 2)
    {
      if (basisDegree == 1)
        return 2;
      else if (basisDegree == 2)
        return 3;
      else if (basisDegree == 3)
        return 4;
    }
    else if (NDIM ==3)
    {
      if (basisDegree == 1)
        return 4;
      else if (basisDegree == 2)
        return 8;
      else if (basisDegree == 3)
        return 12;
      else if (basisDegree == 4)
        return 17;
    }
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
    unsigned nodesPerSide = basisDegree + 1;
    // Number of nodes on a 2-D serendipity element of same degree
    unsigned nodesPerCell2D;
    // Figure out number of nodes per row
    unsigned nodesPerRow2D[nodesPerSide];
    unsigned globalNodesPerRow2D[nodesPerSide];
    if (basisDegree == 1)
    {
      nodesPerCell2D = 4;
      nodesPerRow2D[0] = 2;
      nodesPerRow2D[1] = 2;
    }
    else if (basisDegree == 2)
    {
      nodesPerCell2D = 8;
      nodesPerRow2D[0] = 3;
      nodesPerRow2D[1] = 2;
      nodesPerRow2D[2] = 3;
    }
    else if (basisDegree == 3)
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
      if (basisDegree == 1)
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
      else if (basisDegree == 2)
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
      else if (basisDegree == 3)
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
      if (basisDegree == 1)
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
      else if (basisDegree == 2)
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
      else if (basisDegree == 3)
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
      else if (basisDegree == 4)
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
      if (basisDegree == 1)
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
      else if (basisDegree == 2)
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
      else if (basisDegree == 3)
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
      if (basisDegree == 1)
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
      else if (basisDegree == 2)
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
      else if (basisDegree == 3)
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
      else if (basisDegree == 4)
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
    double coordScales[NDIM];
    for (unsigned i = 0; i < NDIM; i++)
    {
      coordScales[i] = 0.5*dq[i];
    }

    grid.getCentroid(xc);

    // Loop over all node locations on reference element and convert them
    // to appropriate coordinates
    for (unsigned i = 0; i < this->getNumNodes(); i++)
    {
      for (unsigned dim = 0; dim < NDIM; dim++)
      {
        // Node list already has the coordinates of nodes on reference element centered at origin
        nodeCoords(i,dim) = xc[dim] + nodeList(i,dim)*coordScales[dim];
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
    for (unsigned i = 0; i < refMass.rows(); i++)
      for (unsigned j = 0; j < refMass.cols(); j++)
        NjNk(i,j) = refMass(i,j);
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getLowerFaceMassMatrix(unsigned dir,
    Lucee::Matrix<double>& NjNk) const
  {
    for (unsigned i = 0; i < refFaceMassLower[dir].rows(); i++)
      for (unsigned j = 0; j < refFaceMassLower[dir].cols(); j++)
        NjNk(i,j) = refFaceMassLower[dir](i,j);
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getUpperFaceMassMatrix(unsigned dir,
    Lucee::Matrix<double>& NjNk) const
  {
    for (unsigned i = 0; i < refFaceMassUpper[dir].rows(); i++)
      for (unsigned j = 0; j < refFaceMassUpper[dir].cols(); j++)
        NjNk(i,j) = refFaceMassUpper[dir](i,j);
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getStiffnessMatrix(Lucee::Matrix<double>& DNjDNk) const
  {
    for (unsigned i = 0; i < refStiffness.rows(); i++)
      for (unsigned j = 0; j < refStiffness.cols(); j++)
        DNjDNk(i,j) = refStiffness(i,j);
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getGradStiffnessMatrix(
    unsigned dir, Lucee::Matrix<double>& DNjNk) const
  {
    for (unsigned i = 0; i < refGradStiffness[dir].rows(); i++)
      for (unsigned j = 0; j < refGradStiffness[dir].cols(); j++)
        DNjNk(i,j) = refGradStiffness[dir](i,j);
  }

  template <unsigned NDIM>
  unsigned
  SerendipityElement<NDIM>::getNumGaussNodes() const
  {
    unsigned totalNodes = 1;

    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      totalNodes *= numGaussPoints;
    }

    return totalNodes;
  }

  template <unsigned NDIM>
  unsigned
  SerendipityElement<NDIM>::getNumSurfGaussNodes() const
  {
    if (NDIM == 2)
    {
      return numGaussPoints;
    }
    else if (NDIM == 3)
    {
      return numGaussPoints*numGaussPoints;
    }
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
      for (int functionIndex = 0; functionIndex < functionEvaluations.rows(); functionIndex++)
        interpMat(gaussIndex,functionIndex) = functionEvaluations(functionIndex,gaussIndex);
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
      for (int functionIndex = 0; functionIndex < lowerSurfaceEvaluations[dir].cols(); functionIndex++)
        interpMat(gaussIndex,functionIndex) = lowerSurfaceEvaluations[dir](gaussIndex,functionIndex);
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
      for (int functionIndex = 0; functionIndex < upperSurfaceEvaluations[dir].cols(); functionIndex++)
        interpMat(gaussIndex,functionIndex) = upperSurfaceEvaluations[dir](gaussIndex,functionIndex);
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getMomentMatrix(unsigned p, Lucee::Matrix<double>& momMat) const
  {
    if (p > 2)
    {
      // moments higher than 2 not supported for now
      Lucee::Except lce("SerendipityElement::getMomentMatrix: Moment matrix of order ");
      lce << p << " not supported";
      throw lce;
    }
    else
    {
      for (int rowIndex = 0; rowIndex < momMatrix[p].rows(); rowIndex++)
      {
        for (int colIndex = 0; colIndex < momMatrix[p].cols(); colIndex++)
        {
          momMat(rowIndex,colIndex) = momMatrix[p](rowIndex,colIndex);
        }
      }
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
    if (basisDegree == 1)
    {
      data[0] = fldPtr[0];
      data[1] = fldPtr_r[0];
      data[2] = fldPtr_tr[0];
      data[3] = fldPtr_t[0];
    }
    else if (basisDegree == 2)
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
    else if (basisDegree == 3)
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
    else if (basisDegree == 4)
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
      Lucee::Except lce("SerendipityElement::extractFromField: basisDegree ");
      lce << basisDegree << " not supported";
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
 
    // Compute maximum polynomial power we want (no transformations)
    // This will be the maximum power of a single variable in anything
    // evaluated (e.g. product of basis functions)
    maxPower = 2*basisDegree;
    // Populate nodeList according to basisDegree
    nodeList = MatrixXd(this->getNumNodes(),NDIM);
    getNodeList(nodeList);

    std::vector<blitz::Array<double,3> > functionVector;
    computeBasisFunctions(functionVector);

    // Compute gaussian quadrature weights and locations in 1-D
    numGaussPoints = (unsigned)((maxPower+1)/2.0 + 0.5);
    gaussPoints  = std::vector<double>(numGaussPoints);
    gaussWeights = std::vector<double>(numGaussPoints);
    legendre_set(numGaussPoints, &gaussPoints[0], &gaussWeights[0]);
    // Fill out gaussian integration nodes for integration over cell volume
    VectorXd gaussNodeVec(NDIM);
    int currNode = 0;
    double totalWeight;
    
    // Figure out how many gaussian points there are
    int totalVolumeGaussNodes = 1.0;
    int totalSurfaceGaussNodes = 1.0;
    for (unsigned i = 0; i < NDIM; i++)
    {
      totalVolumeGaussNodes = gaussPoints.size()*totalVolumeGaussNodes;
      if (i < NDIM-1)
        totalSurfaceGaussNodes = gaussPoints.size()*totalSurfaceGaussNodes;
    }
    
    gaussNodeList = Eigen::MatrixXd(totalVolumeGaussNodes,NDIM+1);
    functionEvaluations  = Eigen::MatrixXd(this->getNumNodes(),totalVolumeGaussNodes);
/** 3D Array containing derivatives of basis functions (rows) evaluated at gaussian integration locations (cols)
    Correspondance between column and gaussian node set is kept track of in gaussNodeList */
    blitz::Array<double,3> functionDEvaluations(this->getNumNodes(),totalVolumeGaussNodes,NDIM);
    
    gaussNodeListUpperSurf  = std::vector<Eigen::MatrixXd>(NDIM);
    gaussNodeListLowerSurf  = std::vector<Eigen::MatrixXd>(NDIM);
    upperSurfaceEvaluations = std::vector<Eigen::MatrixXd>(NDIM);
    lowerSurfaceEvaluations = std::vector<Eigen::MatrixXd>(NDIM);
    
    if (NDIM == 2)
    {
      // Compute surface quadrature locations
      for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
      {
        // Initialize matrices
        gaussNodeListUpperSurf[dimIndex]  = Eigen::MatrixXd(totalSurfaceGaussNodes,NDIM+1);
        gaussNodeListLowerSurf[dimIndex]  = Eigen::MatrixXd(totalSurfaceGaussNodes,NDIM+1);
        upperSurfaceEvaluations[dimIndex] = Eigen::MatrixXd(totalSurfaceGaussNodes,functionVector.size());
        lowerSurfaceEvaluations[dimIndex] = Eigen::MatrixXd(totalSurfaceGaussNodes,functionVector.size());
        
        // Evaluate all basis functions at all upper and lower surfaces
        for (int nodeIndex = 0; nodeIndex < gaussPoints.size(); nodeIndex++)
        {
          gaussNodeVec(0) = gaussPoints[nodeIndex];
          gaussNodeVec(1) = gaussPoints[nodeIndex];
          gaussNodeVec(dimIndex) = 1;
          totalWeight = gaussWeights[nodeIndex];
          gaussNodeListUpperSurf[dimIndex].row(nodeIndex) << gaussNodeVec(0),gaussNodeVec(1),totalWeight;
          // Evaluate all basis functions at this node location
          for (int functionIndex = 0; functionIndex < functionVector.size(); functionIndex++)
            upperSurfaceEvaluations[dimIndex](nodeIndex, functionIndex) = evalPolynomial(functionVector[functionIndex],gaussNodeVec);
          
          gaussNodeVec(dimIndex) = -1;
          gaussNodeListLowerSurf[dimIndex].row(nodeIndex) << gaussNodeVec(0),gaussNodeVec(1),totalWeight;
          // Evaluate all basis functions at this node location
          for (int functionIndex = 0; functionIndex < functionVector.size(); functionIndex++)
            lowerSurfaceEvaluations[dimIndex](nodeIndex, functionIndex) = evalPolynomial(functionVector[functionIndex],gaussNodeVec);
        }
      }

      // Evaluate all basis functions and their derivatives in every direction at all gaussian volume nodes
      for (unsigned yNodeIndex = 0; yNodeIndex < gaussPoints.size(); yNodeIndex++)
      {
        for (unsigned xNodeIndex = 0; xNodeIndex < gaussPoints.size(); xNodeIndex++)
        {
          gaussNodeVec(0)    = gaussPoints[xNodeIndex];
          gaussNodeVec(1)    = gaussPoints[yNodeIndex];
          totalWeight        = gaussWeights[xNodeIndex]*gaussWeights[yNodeIndex];
          gaussNodeList.row(currNode) << gaussNodeVec(0),gaussNodeVec(1),totalWeight;
          
          for (unsigned functionIndex = 0; functionIndex < functionVector.size(); functionIndex++)
          {
            functionEvaluations(functionIndex,currNode) = evalPolynomial(functionVector[functionIndex],gaussNodeVec);
            for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
            {
              functionDEvaluations((int)functionIndex,currNode,dimIndex) = evalPolynomial(computePolynomialDerivative(functionVector[functionIndex],dimIndex),gaussNodeVec);
            }
          }
          currNode++;
        }
      }
    }
    else if (NDIM == 3)
    {
      // Compute surface quadrature locations
      for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
      {
        // Initialize matrices
        gaussNodeListUpperSurf[dimIndex]  = Eigen::MatrixXd(totalSurfaceGaussNodes,NDIM+1);
        gaussNodeListLowerSurf[dimIndex]  = Eigen::MatrixXd(totalSurfaceGaussNodes,NDIM+1);
        upperSurfaceEvaluations[dimIndex] = Eigen::MatrixXd(totalSurfaceGaussNodes,functionVector.size());
        lowerSurfaceEvaluations[dimIndex] = Eigen::MatrixXd(totalSurfaceGaussNodes,functionVector.size());
        
        // Evaluate all basis functions at all upper and lower surfaces
        currNode = 0;
        for (unsigned dir1NodeIndex = 0; dir1NodeIndex < gaussPoints.size(); dir1NodeIndex++)
        {
          for (unsigned dir2NodeIndex = 0; dir2NodeIndex < gaussPoints.size(); dir2NodeIndex++)
          {
            if (dimIndex == 0)
            {
              gaussNodeVec(0) = 1;
              gaussNodeVec(1) = gaussPoints[dir1NodeIndex];
              gaussNodeVec(2) = gaussPoints[dir2NodeIndex];
              totalWeight = gaussWeights[dir1NodeIndex]*gaussWeights[dir2NodeIndex];
            }
            else if (dimIndex == 1)
            {
              gaussNodeVec(0) = gaussPoints[dir1NodeIndex];
              gaussNodeVec(1) = 1;
              gaussNodeVec(2) = gaussPoints[dir2NodeIndex];
              totalWeight = gaussWeights[dir1NodeIndex]*gaussWeights[dir2NodeIndex];
            }
            else if (dimIndex == 2)
            {
              gaussNodeVec(0) = gaussPoints[dir1NodeIndex];
              gaussNodeVec(1) = gaussPoints[dir2NodeIndex];
              gaussNodeVec(2) = 1;
              totalWeight = gaussWeights[dir1NodeIndex]*gaussWeights[dir2NodeIndex];
            }
          
            gaussNodeListUpperSurf[dimIndex].row(currNode) << gaussNodeVec(0),gaussNodeVec(1),gaussNodeVec(2),totalWeight;
            
            // Evaluate all basis functions at this node location
            for (int functionIndex = 0; functionIndex < functionVector.size(); functionIndex++)
              upperSurfaceEvaluations[dimIndex](currNode, functionIndex) = evalPolynomial(functionVector[functionIndex],gaussNodeVec);
            
            gaussNodeVec(dimIndex) = -1;
            gaussNodeListLowerSurf[dimIndex].row(currNode) << gaussNodeVec(0),gaussNodeVec(1),gaussNodeVec(2),totalWeight;
            // Evaluate all basis functions at this node location
            for (int functionIndex = 0; functionIndex < functionVector.size(); functionIndex++)
              lowerSurfaceEvaluations[dimIndex](currNode, functionIndex) = evalPolynomial(functionVector[functionIndex],gaussNodeVec);

            currNode++;
          }
        }
      }

      // Evaluate all basis functions and their derivatives in every direction at all gaussian volume nodes
      currNode = 0;
      for (unsigned xNodeIndex = 0; xNodeIndex < gaussPoints.size(); xNodeIndex++)
      {
        for (unsigned yNodeIndex = 0; yNodeIndex < gaussPoints.size(); yNodeIndex++)
        {
          for (unsigned zNodeIndex = 0; zNodeIndex < gaussPoints.size(); zNodeIndex++)
          {
            gaussNodeVec(0)    = gaussPoints[xNodeIndex];
            gaussNodeVec(1)    = gaussPoints[yNodeIndex];
            gaussNodeVec(2)    = gaussPoints[zNodeIndex];
            totalWeight        = gaussWeights[xNodeIndex]*gaussWeights[yNodeIndex]*gaussWeights[zNodeIndex];
            gaussNodeList.row(currNode) << gaussNodeVec(0),gaussNodeVec(1),gaussNodeVec(2),totalWeight;
            
            for (unsigned functionIndex = 0; functionIndex < functionVector.size(); functionIndex++)
            {
              functionEvaluations(functionIndex,currNode) = evalPolynomial(functionVector[functionIndex],gaussNodeVec);
              for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
                functionDEvaluations((int)functionIndex,currNode,dimIndex) = evalPolynomial(computePolynomialDerivative(functionVector[functionIndex],dimIndex),gaussNodeVec);
            }
            currNode++;
          }
        }
      }
    }

    // Resize the output matrices we need
    resizeMatrices();
    // Call various functions to populate the matrices

    if (NDIM == 2)
    {
      setupMomentMatrices();
      for(int pIndex = 0; pIndex < 3; pIndex++)
      {
        // Scale into physical space
        for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
          momMatrix[pIndex] *= 0.5*dq[dimIndex];
      }
      
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
    }
    
    computeMass(refMass);
    computeStiffness(functionDEvaluations,refStiffness);

    for (unsigned dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      computeFaceMass(functionVector,dimIndex,refFaceMassLower[dimIndex],refFaceMassUpper[dimIndex]);
      computeGradStiffness(functionDEvaluations,dimIndex,refGradStiffness[dimIndex]);
    }
    
    //printAllMatrices();

    // Scale the matrices computed on reference element into physical space
    for (unsigned dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      refMass   *= 0.5*dq[dimIndex];
      refStiffness *= 0.5*dq[dimIndex];

      // Scale face-mass matrices
      for (unsigned matrixIndex = 0; matrixIndex < NDIM; matrixIndex++)
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
    unsigned generalDim = this->getNumNodes();

    refMass          = Eigen::MatrixXd(generalDim,generalDim);
    refStiffness     = Eigen::MatrixXd(generalDim,generalDim);

    refFaceMassLower = std::vector<Eigen::MatrixXd>(NDIM);
    refFaceMassUpper = std::vector<Eigen::MatrixXd>(NDIM);
    refGradStiffness = std::vector<Eigen::MatrixXd>(NDIM);

    for (unsigned dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      refFaceMassLower[dimIndex] = Eigen::MatrixXd(generalDim,getNumSurfLowerNodes(dimIndex));
      refFaceMassUpper[dimIndex] = Eigen::MatrixXd(generalDim,getNumSurfUpperNodes(dimIndex));
      refGradStiffness[dimIndex] = Eigen::MatrixXd(generalDim,generalDim);
    }
  }
 
  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getNodeList(Eigen::MatrixXd& nodeMatrix)
  {
    if (NDIM == 2)
    {
      if (basisDegree == 1)
      {
        nodeMatrix << -1,-1,
                      1,-1,
                      1,1,
                      -1,1;
      }
      else if (basisDegree == 2)
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
      else if (basisDegree = 3)
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
      if (basisDegree == 1) {
        nodeMatrix << -1,-1,-1,
                     1,-1,-1,
                     1,1,-1,
                     -1,1,-1,
                     -1,-1,1,
                     1,-1,1,
                     1,1,1,
                     -1,1,1;
      }
      else if (basisDegree == 2)
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
      else if (basisDegree == 3)
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
      else if (basisDegree == 4)
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
    // Superlinear degree
    unsigned superDegree;
    unsigned dataIndex = 0;

    if (NDIM == 2)
    {
      for (unsigned yIndex = 0; yIndex <= basisDegree; yIndex++)
      {
        for (unsigned xIndex = 0; xIndex <= basisDegree; xIndex++)
        {
          superDegree = (xIndex>1)*xIndex + (yIndex>1)*yIndex;
          
          if (superDegree < basisDegree + 1)
          {
            basisMatrix.row(dataIndex) = Vector2i(xIndex,yIndex);
            dataIndex = dataIndex + 1;
          }
        }
      }
    }
    else if (NDIM == 3)
    {
      for (unsigned zIndex = 0; zIndex <= basisDegree; zIndex++)
      {
        for (unsigned yIndex = 0; yIndex <= basisDegree; yIndex++)
        {
          for (unsigned xIndex = 0; xIndex <= basisDegree; xIndex++)
          {
            superDegree = (xIndex>1)*xIndex + (yIndex>1)*yIndex + (zIndex>1)*zIndex;
            
            if (superDegree < basisDegree + 1)
            {
              basisMatrix.row(dataIndex) = Vector3i(xIndex,yIndex,zIndex);
              dataIndex = dataIndex + 1;
            }
          }
        }
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::computeBasisFunctions(std::vector<blitz::Array<double,3> >& functionVector)
  {

    // Populate basisList according to basisDegree
    // Matrix to represent basis monomials
    Eigen::MatrixXi basisList(this->getNumNodes(),NDIM);
    setupBasisMatrix(basisList);
    
    // Compute coefficients for basis functions
    MatrixXd coeffMatrix(this->getNumNodes(),this->getNumNodes());
    double coeffProduct;

    // Compute monomial terms evaluated at each nodal coordinate
    for (unsigned nodeIndex = 0; nodeIndex < nodeList.rows(); nodeIndex++)
    {
      for (unsigned basisIndex = 0; basisIndex < basisList.rows(); basisIndex++)
      {
        coeffProduct = 1.0;
        for (unsigned dimIndex = 0; dimIndex < basisList.cols(); dimIndex++)
        {
          coeffProduct = coeffProduct*pow(nodeList(nodeIndex,dimIndex),basisList(basisIndex,dimIndex));
        }
        coeffMatrix(nodeIndex,basisIndex) = coeffProduct;
      }
    }

    // Each column of invCoeffMatrix will be coefficient of basis monomial in basisList
    MatrixXd invCoeffMatrix = coeffMatrix.inverse();

    blitz::Array<double,3> polynomial3DArray;

    if (NDIM == 2)
    {
      blitz::Array<double,3> shapeArray(maxPower,maxPower,1);
      polynomial3DArray.resize(shapeArray.shape());
      polynomial3DArray = 0;
    }
    else if (NDIM == 3)
    {
      blitz::Array<double,3> shapeArray(maxPower,maxPower,maxPower);
      polynomial3DArray.resize(shapeArray.shape());
      polynomial3DArray = 0;
    }

    int xPow,yPow,zPow;

    // Loop over the coefficient matrix (each column represents a basis function)
    for (int functionIndex = 0; functionIndex < invCoeffMatrix.cols(); functionIndex++)
    {
      // Reset the polynomial3DArray
      polynomial3DArray = 0;
      // Loop over the monomial terms that make up the basis function
      for (int monomialIndex = 0; monomialIndex < invCoeffMatrix.rows(); monomialIndex++)
      {
        // Store the basis function in polynomial3DArray by storing its coefficient at
        // a location specified by the monimial term powers
        xPow = basisList(monomialIndex,0);
        yPow = basisList(monomialIndex,1);
        if (NDIM == 2)
          zPow = 0;
        else if (NDIM == 3)
          zPow = basisList(monomialIndex,2);
        polynomial3DArray(xPow,yPow,zPow) = polynomial3DArray(xPow,yPow,zPow) + invCoeffMatrix(monomialIndex,functionIndex);
      }
      // Store the newly represented basis function into functionVector
      functionVector.push_back(polynomial3DArray.copy());
    }
  }

  template <unsigned NDIM>
  double
  SerendipityElement<NDIM>::evalPolynomial(const blitz::Array<double,3>& polyCoeffs, const VectorXd& nodeCoords)
  {
    double totalSum = 0.0;
    double monomialTerm;

    // TODO: Combine NDIM == 2 with NDIM == 3 case by looping over pow(nodeCoords(i),iIndex) multiplications
    if (NDIM == 2)
    {
      for (int xIndex = 0; xIndex < polyCoeffs.extent(0); xIndex++)
      {
        for (int yIndex = 0; yIndex < polyCoeffs.extent(1); yIndex++)
        {
          // Maybe add a check to see if polyCoeffs is 0 so we don't have to compute pow's
          monomialTerm = pow(nodeCoords(0),xIndex)*pow(nodeCoords(1),yIndex);
          totalSum += polyCoeffs(xIndex,yIndex,0)*monomialTerm;
        }
      }
    }
    else if (NDIM == 3)
    {
      for (int xIndex = 0; xIndex < polyCoeffs.extent(0); xIndex++)
      {
        for (int yIndex = 0; yIndex < polyCoeffs.extent(1); yIndex++)
        {
          for (int zIndex = 0; zIndex < polyCoeffs.extent(2); zIndex++)
          {
            // Maybe add a check to see if polyCoeffs is 0 so we don't have to compute pow's
            monomialTerm = pow(nodeCoords(0),xIndex)*pow(nodeCoords(1),yIndex)*pow(nodeCoords(2),zIndex);
            totalSum += polyCoeffs(xIndex,yIndex,zIndex)*monomialTerm;
          }
        }
      }
    }
    return totalSum;
  }
  
  template <unsigned NDIM>
  blitz::Array<double,3>
  SerendipityElement<NDIM>::computePolynomialProduct(const blitz::Array<double,3>& poly1, 
    const blitz::Array<double,3>& poly2)
  {
    int x3Index;
    int y3Index;
    int z3Index;

    blitz::Array<double,3> poly3(poly1.extent(0),poly1.extent(1),poly1.extent(2));
    poly3 = 0;
    
    // Loop over each element of poly1, multiplying it by every element of poly2...
    for (int x1Index = 0; x1Index < poly1.extent(0); x1Index++)
    {
      for (int y1Index = 0; y1Index < poly1.extent(1); y1Index++)
      {
        for (int z1Index = 0; z1Index < poly1.extent(2); z1Index++)
        {
          for (int x2Index = 0; x2Index < poly2.extent(0); x2Index++)
          {
            for (int y2Index = 0; y2Index < poly2.extent(1); y2Index++)
            {
              for (int z2Index = 0; z2Index < poly2.extent(2); z2Index++)
              {
                x3Index = x1Index + x2Index;
                y3Index = y1Index + y2Index;
                z3Index = z1Index + z2Index;
                
                if ((x3Index > poly3.extent(0)) || (y3Index > poly3.extent(1)) || (z3Index > poly3.extent(2)))
                  continue;

                poly3(x3Index,y3Index,z3Index) += poly1(x1Index,y1Index,z1Index)*poly2(x2Index,y2Index,z2Index);
              }
            }
          }
        }
      }
    }

    return poly3;
  }

  template <unsigned NDIM>
  blitz::Array<double,3>
  SerendipityElement<NDIM>::computePolynomialDerivative(const blitz::Array<double,3>& poly, unsigned dir)
  {
    int xResultIndex;
    int yResultIndex;
    int zResultIndex;
    double resultCoefficient;

    blitz::Array<double,3> polyResult(poly.extent(0),poly.extent(1),poly.extent(2));
    polyResult = 0;
    
    // Loop over each element of poly1, multiplying it by every element of poly2...
    for (int xIndex = 0; xIndex < poly.extent(0); xIndex++)
    {
      for (int yIndex = 0; yIndex < poly.extent(1); yIndex++)
      {
        for (int zIndex = 0; zIndex < poly.extent(2); zIndex++)
        {
          if (dir == 0)
          {
            // Compute new location of coefficient
            if (xIndex == 0)
              continue;
            else {
              xResultIndex = xIndex - 1;
              yResultIndex = yIndex;
              zResultIndex = zIndex;
            }
            // Compute and store new coefficient
            polyResult(xResultIndex,yResultIndex,zResultIndex) += xIndex*poly(xIndex,yIndex,zIndex);
          }
          else if (dir == 1)
          {
            // Compute new location of coefficient
            if (yIndex == 0)
              continue;
            else {
              xResultIndex = xIndex;
              yResultIndex = yIndex - 1;
              zResultIndex = zIndex;
            }
            // Compute and store new coefficient
            polyResult(xResultIndex,yResultIndex,zResultIndex) += yIndex*poly(xIndex,yIndex,zIndex);
          }
          else if (dir == 2)
          {
            // Compute new location of coefficient
            if (zIndex == 0)
              continue;
            else {
              xResultIndex = xIndex;
              yResultIndex = yIndex;
              zResultIndex = zIndex-1;
            }
            // Compute and store new coefficient
            polyResult(xResultIndex,yResultIndex,zResultIndex) += zIndex*poly(xIndex,yIndex,zIndex);
          }
        }
      }
    }

    return polyResult;
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::computeMass(Eigen::MatrixXd& resultMatrix)
  {
    blitz::Array<double,2> resultArray(this->getNumNodes(),this->getNumNodes());
    double integrationResult;

    resultMatrix.setZero(resultMatrix.rows(),resultMatrix.cols());
    
    for (unsigned kIndex = 0; kIndex < resultMatrix.rows(); kIndex++)
    {
      for (unsigned mIndex = 0; mIndex < resultMatrix.cols(); mIndex++)
      {
        // Reset integration result
        integrationResult = 0.0;

        // Loop over 3d gauss points to evaluate integral using gaussian quadrature
        for (unsigned gaussIndex = 0; gaussIndex < gaussNodeList.rows(); gaussIndex++)
        {
          // result = weight(node) * f1(node) * f2(node)
          integrationResult += gaussNodeList(gaussIndex,NDIM)*functionEvaluations(kIndex,gaussIndex)*functionEvaluations(mIndex,gaussIndex);
        }
        
        resultMatrix(kIndex,mIndex) = integrationResult;
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::computeFaceMass(const std::vector<blitz::Array<double,3> >& functionVector, unsigned dir,
    Eigen::MatrixXd& lowerResultMatrix, Eigen::MatrixXd& upperResultMatrix)
  {
    std::vector<int> surfLowerNodeNums(lowerResultMatrix.cols(),0);
    std::vector<int> surfUpperNodeNums(upperResultMatrix.cols(),0);

    getSurfLowerNodeNums(dir,surfLowerNodeNums);
    getSurfUpperNodeNums(dir,surfUpperNodeNums);
    
    double integrationResultU;
    double integrationResultL;
    unsigned upperNodeNum;
    unsigned lowerNodeNum;

    for (unsigned kIndex = 0; kIndex < upperResultMatrix.rows(); kIndex++)
    {
      // Note that mIndex will be those of the ones in lower/upper node nums
      for (unsigned mIndex = 0; mIndex < upperResultMatrix.cols(); mIndex++)
      {
        lowerNodeNum = surfLowerNodeNums[mIndex];
        upperNodeNum = surfUpperNodeNums[mIndex];
        integrationResultL = 0.0;
        integrationResultU = 0.0;
        
        for (int testIndex = 0; testIndex < gaussNodeListUpperSurf[dir].rows(); testIndex++)
        {
          integrationResultU += upperSurfaceEvaluations[dir](testIndex,upperNodeNum)*
            upperSurfaceEvaluations[dir](testIndex,kIndex)*gaussNodeListUpperSurf[dir](testIndex,NDIM);
          integrationResultL += lowerSurfaceEvaluations[dir](testIndex,lowerNodeNum)*
            lowerSurfaceEvaluations[dir](testIndex,kIndex)*gaussNodeListLowerSurf[dir](testIndex,NDIM);
        }

        upperResultMatrix(kIndex,mIndex) = integrationResultU;
        lowerResultMatrix(kIndex,mIndex) = integrationResultL;
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::computeStiffness(const blitz::Array<double,3>& functionDerivative,
    Eigen::MatrixXd& resultMatrix)
  {
    double integrationResult;

    resultMatrix.setZero(resultMatrix.rows(),resultMatrix.cols());
    
    for (int kIndex = 0; kIndex < resultMatrix.rows(); kIndex++)
    {
      for (int mIndex = 0; mIndex < resultMatrix.cols(); mIndex++)
      {
        // Reset integration result
        integrationResult = 0.0;
        
        // Loop over 3d gauss points to evaluate integral using gaussian quadrature
        for (int gaussIndex = 0; gaussIndex < gaussNodeList.rows(); gaussIndex++)
        {
          for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
          {
            integrationResult += gaussNodeList(gaussIndex,NDIM)*functionDerivative(kIndex,gaussIndex,dimIndex)*functionDerivative(mIndex,gaussIndex,dimIndex)*4.0/dq2[dimIndex];
          }
        }
        resultMatrix(kIndex,mIndex) = integrationResult;
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::computeGradStiffness(const blitz::Array<double,3>& functionDerivative,
    unsigned dir, Eigen::MatrixXd& resultMatrix)
  {
    double integrationResult;

    resultMatrix.setZero(resultMatrix.rows(),resultMatrix.cols());

    for (unsigned kIndex = 0; kIndex < resultMatrix.rows(); kIndex++)
    {
      for (unsigned mIndex = 0; mIndex < resultMatrix.cols(); mIndex++)
      {
        // Reset integration result
        integrationResult = 0.0;

        // Loop over 3d gauss points to evaluate integral using gaussian quadrature
        for (int gaussIndex = 0; gaussIndex < gaussNodeList.rows(); gaussIndex++)
        {
          integrationResult += gaussNodeList(gaussIndex,NDIM)*functionDerivative((int)kIndex,gaussIndex,(int)dir)*functionEvaluations(mIndex,gaussIndex);
        }

        resultMatrix(kIndex,mIndex) = integrationResult;
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

    // Currently computing moments up to p = 2;
    unsigned momentMax = 2;
    momMatrix = std::vector<Eigen::MatrixXd>(momentMax+1);

    // Initialize Matrices
    for (int momentVal = 0; momentVal < momentMax + 1; momentVal++)
      momMatrix[momentVal] = Eigen::MatrixXd::Zero(nodesPerSide,functionEvaluations.rows());

    // Evaluate Integral y^p*phi(x)psi(x,y)dA
    for (int j = 0; j < nodesPerSide; j++)
    {
      // 1-d basis function in this loop
      int lowerNodeNum = surfLowerNodeNums[j];
      for (int k = 0; k < functionEvaluations.rows(); k++)
      {
        for (int nodeIndex = 0; nodeIndex < gaussNodeList.rows(); nodeIndex++)
        {
          double gaussWeight  = gaussNodeList(nodeIndex,NDIM);
          double yCoord       = gaussNodeList(nodeIndex,1);
          unsigned rollingIndex = nodeIndex % nodesPerSide;

          double yPowerFactor = 1.0;
          for (int momentVal = 0; momentVal <= momentMax; momentVal++)
          {
            if (momentVal != 0)
              yPowerFactor *= yCoord;
            momMatrix[momentVal](j,k) += gaussWeight*yPowerFactor*functionEvaluations(k,nodeIndex)*
              lowerSurfaceEvaluations[1](rollingIndex,lowerNodeNum);
          }
        }
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::printAllMatrices()
  {
    std::cout << "refMass " << std::endl << refMass << std::endl;

    for (unsigned dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      std::cout << "refFaceMass_Lower" << dimIndex << std::endl << refFaceMassLower[dimIndex] << std::endl;
      std::cout << "refFaceMass_Upper" << dimIndex << std::endl << refFaceMassUpper[dimIndex] << std::endl;
    }

    std::cout << "refStiffness " << std::endl << refStiffness << std::endl;

    for (unsigned dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      std::cout << "refGradStiffness_" << dimIndex << std::endl << refGradStiffness[dimIndex] << std::endl;
    }
  }

  
// instantiations
  template class SerendipityElement<1>;
  template class SerendipityElement<2>;
  template class SerendipityElement<3>;
  template class SerendipityElement<4>;
  template class SerendipityElement<5>;
}
