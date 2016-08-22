/**
 * @file	LcSerendipityElement.cpp
 *
 * @brief Serendipity element for DG solvers in several dimensions and orders.
 * Additional code may be required for finite element solves.
 * Should be able to do any dimension and order as long as nodes are specified in a list.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcExcept.h>
#include <LcGridIfc.h>
#include <LcSerendipityElement.h>

// std includes
#include <ctime>

//#define TELL(s) s
#define TELL(s)

namespace Lucee
{
  using namespace Eigen;

// set module name
  template <> const char *SerendipityElement<1>::id = "SerendipityElement";
  template <> const char *SerendipityElement<2>::id = "SerendipityElement";
  template <> const char *SerendipityElement<3>::id = "SerendipityElement";
  template <> const char *SerendipityElement<4>::id = "SerendipityElement";
  template <> const char *SerendipityElement<5>::id = "SerendipityElement";

  static
  double getElapsedTime(clock_t& tmStart, clock_t& tmEnd)
  {
    return (double) (tmEnd-tmStart)/CLOCKS_PER_SEC;
  }

  template <unsigned NDIM>
  SerendipityElement<NDIM>::SerendipityElement()
    : NodalFiniteElementIfc<NDIM>(1), idxr(
        &Lucee::FixedVector<2, unsigned>((unsigned)0)[0], &Lucee::FixedVector<2, int>(1)[0])
  {
    // notice funcky initialization of indexer: this is just so code
    // compiles as indexers do not have default ctors. It gets reset in
    // the readInput() method anyway.
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
    
    if (tbl.hasNumber("num1DGaussPoints"))
      numGaussPoints = tbl.getNumber("num1DGaussPoints");
    else numGaussPoints = (unsigned)((maxPower+1)/2.0 + 0.5);
     
    // Local-to-global mapping valid for polyOrder = 1 in 2d
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->template getGrid<Lucee::StructuredGridBase<NDIM> >();
    Lucee::Region<NDIM, int> grgn = grid.getGlobalRegion();
    
    int shape[NDIM];
    
    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
      shape[dimIndex] = grgn.getShape(dimIndex) + 1;
    Lucee::Region<NDIM, int> lgRgn(shape);
    
    idxr = RowMajorIndexer<NDIM>(lgRgn);
    
    if (NDIM == 1)
    {
      if (polyOrder < 11)
        this->setNumNodes(getSerendipityDimension(polyOrder, NDIM));
      else
      {
        Lucee::Except lce("SerendipityElement: Degree must be 10 or less.");
        lce << " Provided " << polyOrder << " instead";
        throw lce;
      }
      setupMatrices();
    }
    else if (NDIM == 2)
    {
      if (polyOrder < 11)
        this->setNumNodes(getSerendipityDimension(polyOrder, NDIM));
      else
      {
        Lucee::Except lce("SerendipityElement: Degree must be 10 or less.");
        lce << " Provided " << polyOrder << " instead";
        throw lce;
      }
      setupMatrices();
    }
    else if (NDIM == 3 || NDIM == 4 || NDIM == 5)
    {
      if (polyOrder < 5)
        this->setNumNodes(getSerendipityDimension(polyOrder, NDIM));
      else
      {
        Lucee::Except lce("SerendipityElement: Degree must be 4 or less");
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
        ndIds[2] = 3;
      }
      else
      {
        Lucee::Except lce("SerendipityElement::getExclusiveNodeIndices polyOrder not supported.");
        lce << " Provided " << polyOrder;
        throw lce;
      }
    }
    else if (NDIM == 3)
    {
      ndIds.clear();
      if (polyOrder == 1)
      {
        ndIds.resize(1);
        ndIds[0] = 0;
      }
      else
      {
        Lucee::Except lce("SerendipityElement::getExclusiveNodeIndices polyOrder not supported.");
        lce << " Provided " << polyOrder;
        throw lce;
      }
    }
    else
    {
      Lucee::Except lce("SerendipityElement::getExclusiveNodeIndices NDIM must be 1, 2, or 3.");
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
  SerendipityElement<NDIM>::F_func(unsigned nx, unsigned ny, int ix, int iy) const
  {
    return (2*nx+1)*iy + (nx+1)*iy + 2*ix;
  }

  template <unsigned NDIM>
  unsigned
  SerendipityElement<NDIM>::G_func(unsigned nx, unsigned ny, int ix, int iy) const
  {
    return (2*nx+1)*iy + (nx+1)*iy + (2*nx+1) + ix;
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getGlobalIndices(int ix, int iy, std::vector<int>& glob,
    std::vector<int>& loc)
  {
    // The code below basically is to compute the local -> global mapping
    // for exclusively owned nodes. Note that on the right edge, top edge
    // and top-right corner one needs to be careful due to the ownership
    // rules of structured grids. This method should probably be promoted
    // to the top-level LcNodalFiniteElementIfc class and replace the
    // getExclusiveNodeIndices, which returns local node numbers and does
    // not account for the special cases for upper edges and corner.
    
    glob.clear();
    loc.clear();

    const Lucee::StructuredGridBase<2>& grid 
      = this->template getGrid<Lucee::StructuredGridBase<2> >();
    Lucee::Region<2, int> grgn = grid.getGlobalRegion();
    int numX = grgn.getShape(0);
    int numY = grgn.getShape(1);

    if ((ix<numX) && (iy<numY))
    { // nodes inside grid proper
      glob.resize(3);
      glob[0] = F_func(numX, numY, ix, iy); // node 0
      glob[1] = F_func(numX, numY, ix, iy) + 1; // node 2
      glob[2] = G_func(numX, numY, ix, iy); // node 4

      loc.resize(3);
      loc[0] = 0;
      loc[1] = 1;
      loc[2] = 2;
    }
    else if ((ix==numX) && (iy<numY))
    { // right edge nodes
      glob.resize(2);
      glob[0] = F_func(numX, numY, ix, iy); // node 0
      glob[1] = G_func(numX, numY, ix, iy); // node 4

      loc.resize(2);
      loc[0] = 0;
      loc[1] = 2;
    }
    else if ((ix<numX) && (iy==numY))
    { // top edge nodes
      glob.resize(2);
      glob[0] = F_func(numX, numY, ix, iy); // node 0
      glob[1] = F_func(numX, numY, ix, iy) + 1; // node 1

      loc.resize(2);
      loc[0] = 0;
      loc[1] = 1;
    }
    else if ((ix==numX) && (iy==numY))
    { // top right corner
      glob.resize(1);
      glob[0] = F_func(numX, numY, ix, iy); // node 0

      loc.resize(1);
      loc[0] = 0;
    }
  }

  template <unsigned NDIM>
  unsigned
  SerendipityElement<NDIM>::getNumGlobalNodes() const
  {
    if (NDIM > 3)
      throw Lucee::Except("SerendipityElement::getNumGlobalNodes: Not implemented for NDIM > 2!");
    // get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->template getGrid<Lucee::StructuredGridBase<NDIM> >();
    Lucee::Region<NDIM, int> gridRgn = grid.getGlobalRegion();
    int numGlobalNodes = 0;

    if (NDIM == 1)
    {
      numGlobalNodes = 1 + polyOrder*gridRgn.getShape(0);
    }
    else if (NDIM == 2)
    {
      int nx = gridRgn.getShape(0);
      int ny = gridRgn.getShape(1);

      if (polyOrder == 1)
        numGlobalNodes = (nx+1)*(ny+1);
      else if (polyOrder ==2)
        // there are 3 owned nodes per cell, 2*nx and 2*ny edge nodes along
        // the top and right edges and the 1 accounts for the node on the
        // top-right corner.
        numGlobalNodes = 3*nx*ny + 2*nx + 2*ny + 1;
    }
    else if (NDIM == 3)
    {
      int nx = gridRgn.getShape(0);
      int ny = gridRgn.getShape(1);
      int nz = gridRgn.getShape(2);

      if (polyOrder == 1)
        numGlobalNodes = (nx+1)*(ny+1)*(nz+1);
      else
      {
        Lucee::Except lce("SerendipityElement::getExclusiveNodeIndices polyOrder not supported.");
        lce << " Provided " << polyOrder;
        throw lce;
      }
    }

    return numGlobalNodes;
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getLocalToGlobal(std::vector<int>& lgMap) const
  {
    if (NDIM == 1)
    {
      // Loop over exclusive nodes
      for (int nodeIndex = 0; nodeIndex < polyOrder; nodeIndex++)
        lgMap[nodeIndex] = this->currIdx[0]*polyOrder + nodeIndex;
    }
    else if (NDIM == 2)
    {
      int ix = this->currIdx[0];
      int iy = this->currIdx[1];

      if (polyOrder == 1)
      {
        lgMap[0] = idxr.getIndex(ix, iy); // 0
        lgMap[1] = idxr.getIndex(ix+1, iy); // 1
        lgMap[2] = idxr.getIndex(ix, iy+1); // 2
        lgMap[3] = idxr.getIndex(ix+1, iy+1); // 3
      }
      else if (polyOrder == 2)
      {
        const Lucee::StructuredGridBase<2>& grid 
        = this->template getGrid<Lucee::StructuredGridBase<2> >();
        Lucee::Region<2, int> grgn = grid.getGlobalRegion();
        int numX = grgn.getShape(0);
        int numY = grgn.getShape(1);

        lgMap[0] = F_func(numX, numY, ix, iy); // 0
        lgMap[1] = F_func(numX, numY, ix, iy) + 1; // 1
        lgMap[2] = F_func(numX, numY, ix, iy) + 2; // 2
        lgMap[3] = G_func(numX, numY, ix, iy); // 3
        lgMap[4] = G_func(numX, numY, ix, iy) + 1; // 4
        lgMap[5] = F_func(numX, numY, ix, iy+1); // 5
        lgMap[6] = F_func(numX, numY, ix, iy+1) + 1; // 6
        lgMap[7] = F_func(numX, numY, ix, iy+1) + 2; // 7
      }
      else throw Lucee::Except("SerendipityElement::getLocalToGlobal: polyOrder not implemented!");
    }
    else if (NDIM == 3)
    {
      int ix = this->currIdx[0];
      int iy = this->currIdx[1];
      int iz = this->currIdx[2];

      if (polyOrder == 1)
      {
        lgMap[0] = idxr.getIndex(ix, iy, iz);
        lgMap[1] = idxr.getIndex(ix+1, iy, iz);
        lgMap[2] = idxr.getIndex(ix, iy+1, iz);
        lgMap[3] = idxr.getIndex(ix+1, iy+1, iz);
        // Same as previous four mappings but with iz+1 instead of iz
        lgMap[4] = idxr.getIndex(ix, iy, iz+1);
        lgMap[5] = idxr.getIndex(ix+1, iy, iz+1);
        lgMap[6] = idxr.getIndex(ix, iy+1, iz+1);
        lgMap[7] = idxr.getIndex(ix+1, iy+1, iz+1);
      }
      else throw Lucee::Except("SerendipityElement::getLocalToGlobal: polyOrder not implemented!");
    }
    else throw Lucee::Except("SerendipityElement::getLocalToGlobal: Dimension not implemented!");

  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getSurfLowerLocalToGlobal(unsigned dir,
    std::vector<int>& lgMap) const
  {
    if (NDIM == 1)
      lgMap[0] = this->currIdx[0]*polyOrder;
    else if (NDIM == 2)
    {
      int ix = this->currIdx[0];
      int iy = this->currIdx[1];

      if (polyOrder == 1)
      {
        if (dir == 0)
        {
          lgMap[0] = idxr.getIndex(ix, iy); // 0
          lgMap[1] = idxr.getIndex(ix, iy+1); // 2
        }
        else if (dir == 1)
        {
          lgMap[0] = idxr.getIndex(ix, iy); // 0
          lgMap[1] = idxr.getIndex(ix+1, iy); // 1
        }
      }
      else if (polyOrder == 2)
      {
        const Lucee::StructuredGridBase<2>& grid 
        = this->template getGrid<Lucee::StructuredGridBase<2> >();
        Lucee::Region<2, int> grgn = grid.getGlobalRegion();
        int numX = grgn.getShape(0);
        int numY = grgn.getShape(1);

        if (dir == 0)
        {
          lgMap[0] = F_func(numX, numY, ix, iy); // 0
          lgMap[1] = G_func(numX, numY, ix, iy); // 3
          lgMap[2] = F_func(numX, numY, ix, iy+1); // 5
        }
        else if (dir == 1)
        {
          lgMap[0] = F_func(numX, numY, ix, iy); // 0
          lgMap[1] = F_func(numX, numY, ix, iy) + 1; // 1
          lgMap[2] = F_func(numX, numY, ix, iy) + 2; // 2
        }
      }
      else throw Lucee::Except("SerendipityElement::getSurfLowerLocalToGlobal: polyOrder not implemented!");
    }
    else if (NDIM == 3)
    {
      int ix = this->currIdx[0];
      int iy = this->currIdx[1];
      int iz = this->currIdx[2];

      if(polyOrder == 1)
      {
        if (dir == 0)
        {
          lgMap[0] = idxr.getIndex(ix, iy, iz); // 0
          lgMap[1] = idxr.getIndex(ix, iy+1, iz); // 2
          lgMap[2] = idxr.getIndex(ix, iy, iz+1); // 4
          lgMap[3] = idxr.getIndex(ix, iy+1, iz+1); // 6
        }
        else if (dir == 1)
        {
          lgMap[0] = idxr.getIndex(ix, iy, iz); // 0
          lgMap[1] = idxr.getIndex(ix+1, iy, iz); // 1
          lgMap[2] = idxr.getIndex(ix, iy, iz+1); // 4
          lgMap[3] = idxr.getIndex(ix+1, iy, iz+1); // 5
        }
        else if (dir == 2)
        {
          lgMap[0] = idxr.getIndex(ix, iy, iz); // 0
          lgMap[1] = idxr.getIndex(ix+1, iy, iz); // 1
          lgMap[2] = idxr.getIndex(ix, iy+1, iz); // 2
          lgMap[3] = idxr.getIndex(ix+1, iy+1, iz); // 3
        }
      }
      else throw Lucee::Except("SerendipityElement::getSurfLowerLocalToGlobal: polyOrder not implemented!");
    }
    else throw Lucee::Except("SerendipityElement::getSurfLowerLocalToGlobal: Dimension not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getSurfUpperLocalToGlobal(unsigned dir,
    std::vector<int>& lgMap) const
  {
    if (NDIM == 1)
      lgMap[0] = (this->currIdx[0]+1)*polyOrder;
    else if (NDIM == 2)
    {
      int ix = this->currIdx[0];
      int iy = this->currIdx[1];

      if (polyOrder == 1)
      {
        if (dir == 0)
        {
          lgMap[0] = idxr.getIndex(ix+1, iy); // 1
          lgMap[1] = idxr.getIndex(ix+1, iy+1); // 3
        }
        else if (dir == 1)
        {
          lgMap[0] = idxr.getIndex(ix, iy+1); // 2
          lgMap[1] = idxr.getIndex(ix+1, iy+1); // 3
        }
      }
      else if (polyOrder == 2)
      {
        const Lucee::StructuredGridBase<2>& grid 
        = this->template getGrid<Lucee::StructuredGridBase<2> >();
        Lucee::Region<2, int> grgn = grid.getGlobalRegion();
        int numX = grgn.getShape(0);
        int numY = grgn.getShape(1);

        if (dir == 0)
        {
          lgMap[0] = F_func(numX, numY, ix, iy) + 2; // 2
          lgMap[1] = G_func(numX, numY, ix, iy) + 1; // 4
          lgMap[2] = F_func(numX, numY, ix, iy+1) + 2; // 7
        }
        else if (dir == 1)
        {
          lgMap[0] = F_func(numX, numY, ix, iy+1); // 5
          lgMap[1] = F_func(numX, numY, ix, iy+1) + 1; // 6
          lgMap[2] = F_func(numX, numY, ix, iy+1) + 2; // 7
        }
      }
      else throw Lucee::Except("SerendipityElement::getSurfUpperLocalToGlobal: polyOrder not implemented!");
    }
    else if (NDIM == 3)
    {
      int ix = this->currIdx[0];
      int iy = this->currIdx[1];
      int iz = this->currIdx[2];

      if(polyOrder == 1)
      {
        if (dir == 0)
        {
          lgMap[0] = idxr.getIndex(ix+1, iy, iz); // 1
          lgMap[1] = idxr.getIndex(ix+1, iy+1, iz); // 3
          lgMap[2] = idxr.getIndex(ix+1, iy, iz+1); // 5
          lgMap[3] = idxr.getIndex(ix+1, iy+1, iz+1); // 7
        }
        else if (dir == 1)
        {
          lgMap[0] = idxr.getIndex(ix, iy+1, iz); // 2
          lgMap[1] = idxr.getIndex(ix+1, iy+1, iz); // 3
          lgMap[2] = idxr.getIndex(ix, iy+1, iz+1); // 6
          lgMap[3] = idxr.getIndex(ix+1, iy+1, iz+1); // 7
        }
        else if (dir == 2)
        {
          lgMap[0] = idxr.getIndex(ix, iy, iz+1); // 4
          lgMap[1] = idxr.getIndex(ix+1, iy, iz+1); // 5
          lgMap[2] = idxr.getIndex(ix, iy+1, iz+1); // 6
          lgMap[3] = idxr.getIndex(ix+1, iy+1, iz+1); // 7
        }
      }
      else throw Lucee::Except("SerendipityElement::getSurfLowerLocalToGlobal: polyOrder not implemented!");
    }
    else throw Lucee::Except("SerendipityElement::getSurfUpperLocalToGlobal: Dimension not implemented!");
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
    if (NDIM == 2)
    {
      if (polyOrder == 1)
      {
        unsigned nn = this->getNumNodes();
        for (int k = 0; k < nn; k++)
          w[k] = 0.5*dq[0];
      }
      else if (polyOrder == 2)
      {
        w[0] = 0.5*dq[0]*(-1.0)/3.0;
        w[1] = 0.5*dq[0]*4.0/3.0;
        w[2] = 0.5*dq[0]*(-1.0)/3.0;
      }
      else if (polyOrder == 3)
      {
        w[0] = 0.5*dq[0]*1.0/4.0;
        w[1] = 0.5*dq[0]*3.0/4.0;
        w[2] = 0.5*dq[0]*3.0/4.0;
        w[3] = 0.5*dq[0]*1.0/4.0;
      }
      else if (polyOrder == 4)
      {
        w[0] = 0.5*dq[0]*7.0/45.0;
        w[1] = 0.5*dq[0]*32.0/45.0;
        w[2] = 0.5*dq[0]*4.0/15.0;
        w[3] = 0.5*dq[0]*32.0/45.0;
        w[4] = 0.5*dq[0]*7.0/45.0;
      }
    }
    else if (NDIM == 2)
    {
      if (polyOrder == 1)
      {
        unsigned nn = this->getNumNodes();
        for (int k = 0; k < nn; k++)
          w[k] = 0.5*dq[0]*0.5*dq[1];
      }
      else if (polyOrder == 2)
      {
        // Lobatto weights
        w[0] = 0.5*dq[0]*0.5*dq[1]*(-1.0)/3.0;
        w[1] = 0.5*dq[0]*0.5*dq[1]*4.0/3.0;
        w[2] = 0.5*dq[0]*0.5*dq[1]*(-1.0)/3.0;
        w[3] = 0.5*dq[0]*0.5*dq[1]*4.0/3.0;
        w[4] = 0.5*dq[0]*0.5*dq[1]*4.0/3.0;
        w[5] = 0.5*dq[0]*0.5*dq[1]*(-1.0)/3.0;
        w[6] = 0.5*dq[0]*0.5*dq[1]*4.0/3.0;
        w[7] = 0.5*dq[0]*0.5*dq[1]*(-1.0)/3.0;
      }
      else if (polyOrder == 3)
      {
        // Lobatto weights
        w[0] = 0.5*dq[0]*0.5*dq[1]*(-1.0)/2.0;
        w[1] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[2] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[3] = 0.5*dq[0]*0.5*dq[1]*(-1.0)/2.0;
        w[4] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[5] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[6] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[7] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[8] = 0.5*dq[0]*0.5*dq[1]*(-1.0)/2.0;
        w[9] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[10] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[11] = 0.5*dq[0]*0.5*dq[1]*(-1.0)/2.0;
      }
      else if (polyOrder == 4)
      {
        // Lobatto weights
        w[0] = 0.5*dq[0]*0.5*dq[1]*(-11.0)/45.0;
        w[1] = 0.5*dq[0]*0.5*dq[1]*32.0/45.0;
        w[2] = 0.5*dq[0]*0.5*dq[1]*(-28.0)/45.0;
        w[3] = 0.5*dq[0]*0.5*dq[1]*32.0/45.0;
        w[4] = 0.5*dq[0]*0.5*dq[1]*(-11.0)/45.0;
        w[5] = 0.5*dq[0]*0.5*dq[1]*32.0/45.0;
        w[6] = 0.5*dq[0]*0.5*dq[1]*32.0/45.0;
        w[7] = 0.5*dq[0]*0.5*dq[1]*(-28.0)/45.0;
        w[8] = 0.5*dq[0]*0.5*dq[1]*16.0/9.0;
        w[9] = 0.5*dq[0]*0.5*dq[1]*(-28.0)/45.0;
        w[10] = 0.5*dq[0]*0.5*dq[1]*32.0/45.0;
        w[11] = 0.5*dq[0]*0.5*dq[1]*32.0/45.0;
        w[12] = 0.5*dq[0]*0.5*dq[1]*(-11.0)/45.0;
        w[13] = 0.5*dq[0]*0.5*dq[1]*32.0/45.0;
        w[14] = 0.5*dq[0]*0.5*dq[1]*(-28.0)/45.0;
        w[15] = 0.5*dq[0]*0.5*dq[1]*32.0/45.0;
        w[16] = 0.5*dq[0]*0.5*dq[1]*(-11.0)/45.0;
      }
    }
    else if (NDIM == 3)
    {
      if (polyOrder == 1)
      {
        unsigned nn = this->getNumNodes();
        for (int k = 0; k < nn; k++)
          w[k] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2];
      }
      else if (polyOrder == 2)
      {
        // Lobatto weights
        w[0] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-1.0);
        w[1] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*4.0/3.0;
        w[2] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-1.0);
        w[3] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*4.0/3.0;
        w[4] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*4.0/3.0;
        w[5] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-1.0);
        w[6] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*4.0/3.0;
        w[7] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-1.0);
        w[8] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*4.0/3.0;
        w[9] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*4.0/3.0;
        w[10] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*4.0/3.0;
        w[11] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*4.0/3.0;
        w[12] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-1.0);
        w[13] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*4.0/3.0;
        w[14] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-1.0);
        w[15] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*4.0/3.0;
        w[16] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*4.0/3.0;
        w[17] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-1.0);
        w[18] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*4.0/3.0;
        w[19] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-1.0);
      }
      else if (polyOrder == 3)
      {
        // Lobatto weights
        w[0] = 0.5*dq[0]*0.5*dq[1]*(-5.0)/4.0;
        w[1] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[2] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[3] = 0.5*dq[0]*0.5*dq[1]*(-5.0)/4.0;
        w[4] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[5] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[6] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[7] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[8] = 0.5*dq[0]*0.5*dq[1]*(-5.0)/4.0;
        w[9] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[10] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[11] = 0.5*dq[0]*0.5*dq[1]*(-5.0)/4.0;
        w[12] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[13] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[14] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[15] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[16] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[17] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[18] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[19] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[20] = 0.5*dq[0]*0.5*dq[1]*(-5.0)/4.0;
        w[21] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[22] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[23] = 0.5*dq[0]*0.5*dq[1]*(-5.0)/4.0;
        w[24] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[25] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[26] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[27] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[28] = 0.5*dq[0]*0.5*dq[1]*(-5.0)/4.0;
        w[29] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[30] = 0.5*dq[0]*0.5*dq[1]*3.0/4.0;
        w[31] = 0.5*dq[0]*0.5*dq[1]*(-5.0)/4.0;

      }
      else if (polyOrder == 4)
      {
        w[0] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-1.0)/5.0;
        w[1] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*32.0/45.0;
        w[2] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-68.0)/45.0;
        w[3] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*32.0/45.0;
        w[4] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-1.0)/5.0;
        w[5] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*32.0/45.0;
        w[6] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*32.0/45.0;
        w[7] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-68.0)/45.0;
        w[8] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*16.0/9.0;
        w[9] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-68.0)/45.0;
        w[10] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*32.0/45.0;
        w[11] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*32.0/45.0;
        w[12] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-1.0)/5.0;
        w[13] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*32.0/45.0;
        w[14] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-68.0)/45.0;
        w[15] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*32.0/45.0;
        w[16] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-1.0)/5.0;
        w[17] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*32.0/45.0;
        w[18] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*32.0/45.0;
        w[19] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*32.0/45.0;
        w[20] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*32.0/45.0;
        w[21] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-68.0)/45.0;
        w[22] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*16.0/9.0;
        w[23] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-68.0)/45.0;
        w[24] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*16.0/9.0;
        w[25] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*16.0/9.0;
        w[26] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-68.0)/45.0;
        w[27] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*16.0/9.0;
        w[28] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-68.0)/45.0;
        w[29] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*32.0/45.0;
        w[30] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*32.0/45.0;
        w[31] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*32.0/45.0;
        w[32] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*32.0/45.0;
        w[33] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-1.0)/5.0;
        w[34] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*32.0/45.0;
        w[35] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-68.0)/45.0;
        w[36] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*32.0/45.0;
        w[37] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-1.0)/5.0;
        w[38] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*32.0/45.0;
        w[39] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*32.0/45.0;
        w[40] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-68.0)/45.0;
        w[41] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*16.0/9.0;
        w[42] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-68.0)/45.0;
        w[43] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*32.0/45.0;
        w[44] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*32.0/45.0;
        w[45] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-1.0)/5.0;
        w[46] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*32.0/45.0;
        w[47] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-68.0)/45.0;
        w[48] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*32.0/45.0;
        w[49] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*(-1.0)/5.0;
      }
    }
    else if (NDIM == 4)
    {
      if (polyOrder == 2)
      {
        unsigned nn = this->getNumNodes();
        for (int k = 0; k < nn; k++)
          w[k] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3];
      }
      else if (polyOrder == 2)
      {
        // Lobatto weights
        w[0] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*(-5.0)/3.0;
        w[1] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[2] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*(-5.0)/3.0;
        w[3] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[4] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[5] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*(-5.0)/3.0;
        w[6] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[7] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*(-5.0)/3.0;
        w[8] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[9] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[10] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[11] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[12] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*(-5.0)/3.0;
        w[13] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[14] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*(-5.0)/3.0;
        w[15] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[16] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[17] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*(-5.0)/3.0;
        w[18] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[19] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*(-5.0)/3.0;
        w[20] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[21] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[22] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[23] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[24] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[25] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[26] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[27] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[28] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*(-5.0)/3.0;
        w[29] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[30] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*(-5.0)/3.0;
        w[31] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[32] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[33] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*(-5.0)/3.0;
        w[34] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[35] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*(-5.0)/3.0;
        w[36] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[37] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[38] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[39] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[40] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*(-5.0)/3.0;
        w[41] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[42] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*(-5.0)/3.0;
        w[43] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[44] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[45] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*(-5.0)/3.0;
        w[46] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*4.0/3.0;
        w[47] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*(-5.0)/3.0;
      }
    }
    else if (NDIM == 5)
    {
      if (polyOrder == 1)
      {
        unsigned nn = this->getNumNodes();
        for (int k = 0; k < nn; k++)
          w[k] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4];
      }
      else if (polyOrder == 2)
      {
        // Lobatto weights
        w[0] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[1] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[2] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[3] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[4] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[5] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[6] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[7] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[8] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[9] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[10] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[11] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[12] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[13] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[14] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[15] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[16] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[17] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[18] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[19] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[20] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[21] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[22] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[23] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[24] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[25] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[26] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[27] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[28] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[29] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[30] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[31] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[32] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[33] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[34] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[35] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[36] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[37] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[38] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[39] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[40] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[41] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[42] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[43] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[44] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[45] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[46] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[47] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[48] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[49] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[50] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[51] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[52] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[53] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[54] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[55] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[56] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[57] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[58] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[59] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[60] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[61] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[62] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[63] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[64] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[65] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[66] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[67] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[68] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[69] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[70] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[71] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[72] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[73] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[74] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[75] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[76] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[77] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[78] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[79] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[80] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[81] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[82] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[83] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[84] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[85] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[86] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[87] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[88] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[89] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[90] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[91] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[92] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[93] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[94] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[95] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[96] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[97] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[98] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[99] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[100] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[101] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[102] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[103] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[104] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[105] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[106] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[107] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[108] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[109] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
        w[110] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*4.0/3.0;
        w[111] = 0.5*dq[0]*0.5*dq[1]*0.5*dq[2]*0.5*dq[3]*0.5*dq[4]*(-7.0)/3.0;
      }
    }
    else throw Lucee::Except("SerendipityElement::getWeights: Dimension not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getSurfUpperWeights(unsigned dir, std::vector<double>& w)
  {
    if (NDIM == 2)
    {
      if (polyOrder == 1)
      {
        unsigned nn = this->getNumSurfUpperNodes(dir);
        if (dir == 0)
        {
          for (int k = 0; k < nn; k++) 
            w[k] = 0.5*dq[1];
        }
        else if (dir == 1)
        {
          for (int k = 0; k < nn; k++) 
            w[k] = 0.5*dq[0];
        }
      }
      else if (polyOrder == 2)
      {
        if (dir == 0)
        {
          w[0] = 0.5*dq[1]/3.0;
          w[1] = 0.5*4*dq[1]/3.0;
          w[2] = 0.5*dq[1]/3.0;
        }
        else if (dir == 1)
        {
          w[0] = 0.5*dq[0]/3.0;
          w[1] = 0.5*4*dq[0]/3.0;
          w[2] = 0.5*dq[0]/3.0;
        }
      }
      else throw Lucee::Except("SerendipityElement::getSurfUpperWeights: polyOrder not implemented!");
    }
    else if (NDIM == 3)
    {
      if (polyOrder == 1)
      {
        unsigned nn = this->getNumSurfUpperNodes(dir);
        if (dir == 0)
        {
          for (int k = 0; k < nn; k++) 
            w[k] = 0.5*dq[1]*0.5*dq[2];
        }
        else if (dir == 1)
        {
          for (int k = 0; k < nn; k++) 
            w[k] = 0.5*dq[0]*0.5*dq[2];
        }
        else if (dir == 2)
        {
          for (int k = 0; k < nn; k++) 
            w[k] = 0.5*dq[0]*0.5*dq[1];
        }
      }
    }
    else throw Lucee::Except("SerendipityElement::getSurfUpperWeights: Dimension not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getSurfLowerWeights(unsigned dir, std::vector<double>& w)
  {
    if (NDIM == 2)
    {
      if (polyOrder == 1)
      {
        unsigned nn = this->getNumSurfUpperNodes(dir);
        if (dir == 0)
        {
          for (int k = 0; k < nn; k++) 
            w[k] = 0.5*dq[1];
        }
        else if (dir == 1)
        {
          for (int k = 0; k < nn; k++) 
            w[k] = 0.5*dq[0];
        }
      }
      else if (polyOrder == 2)
      {
        if (dir == 0)
        {
          w[0] = 0.5*dq[1]/3.0;
          w[1] = 0.5*4*dq[1]/3.0;
          w[2] = 0.5*dq[1]/3.0;
        }
        else if (dir == 1)
        {
          w[0] = 0.5*dq[0]/3.0;
          w[1] = 0.5*4*dq[0]/3.0;
          w[2] = 0.5*dq[0]/3.0;
        }
      }
      else throw Lucee::Except("SerendipityElement::getSurfLowerWeights: polyOrder not implemented!");
    }
    else if (NDIM == 3)
    {
      if (polyOrder == 1)
      {
        unsigned nn = this->getNumSurfUpperNodes(dir);
        if (dir == 0)
        {
          for (int k = 0; k < nn; k++) 
            w[k] = 0.5*dq[1]*0.5*dq[2];
        }
        else if (dir == 1)
        {
          for (int k = 0; k < nn; k++) 
            w[k] = 0.5*dq[0]*0.5*dq[2];
        }
        else if (dir == 2)
        {
          for (int k = 0; k < nn; k++) 
            w[k] = 0.5*dq[0]*0.5*dq[1];
        }
      }
    }
    else throw Lucee::Except("SerendipityElement::getSurfLowerWeights: Dimension not implemented!");
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
  SerendipityElement<NDIM>::getPerpStiffnessMatrix(Lucee::Matrix<double>& DNjDNk) const
  {
    for (int i = 0; i < refPerpStiffness.rows(); i++)
      for (int j = 0; j < refPerpStiffness.cols(); j++)
        DNjDNk(i, j) = refPerpStiffness(i, j);
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
    if (NDIM != 2 && NDIM != 3)
      Lucee::Except("SerendipityElement::getDiffusionMatrices: Only implemented for 2D and 3D");

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
    if (NDIM != 2 && NDIM != 3)
      Lucee::Except("SerendipityElement::getLowerReflectingBcMapping: Only implemented for 2D");

    if (NDIM == 2)
    {
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
          nodeMap[3] = 4;
          nodeMap[4] = 3;
          nodeMap[5] = 7;
          nodeMap[6] = 6;
          nodeMap[7] = 5;
        }
      }
      else if (dir == 1)
      {
        if (polyOrder == 1)
        {
          nodeMap[0] = 2;
          nodeMap[1] = 3;
          nodeMap[2] = 0;
          nodeMap[3] = 1;
        }
        else if (polyOrder == 2)
        {
          nodeMap[0] = 5;
          nodeMap[1] = 6;
          nodeMap[2] = 7;
          nodeMap[3] = 3;
          nodeMap[4] = 4;
          nodeMap[5] = 0;
          nodeMap[6] = 1;
          nodeMap[7] = 2;
        }
      }
    }
    else if (NDIM == 3)
    {
      if (dir == 0)
      {
        if (polyOrder == 1)
        {
          nodeMap[0] = 1;
          nodeMap[1] = 0;
          nodeMap[2] = 3;
          nodeMap[3] = 2;
          nodeMap[4] = 5;
          nodeMap[5] = 4;
          nodeMap[6] = 7;
          nodeMap[7] = 6;
        }
        else if (polyOrder == 2)
        {
          nodeMap[0] = 2;
          nodeMap[1] = 1;
          nodeMap[2] = 0;
          nodeMap[3] = 4;
          nodeMap[4] = 3;
          nodeMap[5] = 7;
          nodeMap[6] = 6;
          nodeMap[7] = 5;
          nodeMap[8] = 9;
          nodeMap[9] = 8;
          nodeMap[10] = 11;
          nodeMap[11] = 10;
          nodeMap[12] = 14;
          nodeMap[13] = 13;
          nodeMap[14] = 12;
          nodeMap[15] = 16;
          nodeMap[16] = 15;
          nodeMap[17] = 19;
          nodeMap[18] = 18;
          nodeMap[19] = 17;
        }
      }
      else if (dir == 1)
      {
        if (polyOrder == 1)
        {
          nodeMap[0] = 2;
          nodeMap[1] = 3;
          nodeMap[2] = 0;
          nodeMap[3] = 1;
          nodeMap[4] = 6;
          nodeMap[5] = 7;
          nodeMap[6] = 4;
          nodeMap[7] = 5;
        }
        else if (polyOrder == 2)
        {
          nodeMap[0] = 5;
          nodeMap[1] = 6;
          nodeMap[2] = 7;
          nodeMap[3] = 3;
          nodeMap[4] = 4;
          nodeMap[5] = 0;
          nodeMap[6] = 1;
          nodeMap[7] = 2;
          nodeMap[8] = 10;
          nodeMap[9] = 11;
          nodeMap[10] = 8;
          nodeMap[11] = 9;
          nodeMap[12] = 17;
          nodeMap[13] = 18;
          nodeMap[14] = 19;
          nodeMap[15] = 15;
          nodeMap[16] = 16;
          nodeMap[17] = 12;
          nodeMap[18] = 13;
          nodeMap[19] = 14;
        }
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getUpperReflectingBcMapping(unsigned dir,
        std::vector<unsigned>& nodeMap) const
  {
    /*
    if (polyOrder > 2)
      Lucee::Except("SerendipityElement::getUpperReflectingBcMapping: Not implemented for higher than quadratic!");
    if (NDIM != 2 && NDIM != 3)
      Lucee::Except("SerendipityElement::getUpperReflectingBcMapping: Only implemented for 2D and 3D");
      */

    for (int srcNodeIndex = 0; srcNodeIndex < nodeList.rows(); srcNodeIndex++)
    {
      for (int refNodeIndex = 0; refNodeIndex < nodeList.rows(); refNodeIndex++)
      {
        if (isReflectionNode(srcNodeIndex, refNodeIndex, dir) == true)
        {
          nodeMap[srcNodeIndex] = refNodeIndex;
          break;
        }
      }
    }
/*
    if (NDIM == 2)
    {
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
          nodeMap[3] = 4;
          nodeMap[4] = 3;
          nodeMap[5] = 7;
          nodeMap[6] = 6;
          nodeMap[7] = 5;
        }
      }
      else if (dir == 1)
      {
        if (polyOrder == 1)
        {
          nodeMap[0] = 2;
          nodeMap[1] = 3;
          nodeMap[2] = 0;
          nodeMap[3] = 1;
        }
        else if (polyOrder == 2)
        {
          nodeMap[0] = 5;
          nodeMap[1] = 6;
          nodeMap[2] = 7;
          nodeMap[3] = 3;
          nodeMap[4] = 4;
          nodeMap[5] = 0;
          nodeMap[6] = 1;
          nodeMap[7] = 2;
        }
      }
    }
    else if (NDIM == 3)
    {
      if (dir == 0)
      {
        if (polyOrder == 1)
        {
          nodeMap[0] = 1;
          nodeMap[1] = 0;
          nodeMap[2] = 3;
          nodeMap[3] = 2;
          nodeMap[4] = 5;
          nodeMap[5] = 4;
          nodeMap[6] = 7;
          nodeMap[7] = 6;
        }
        else if (polyOrder == 2)
        {
          nodeMap[0] = 2;
          nodeMap[1] = 1;
          nodeMap[2] = 0;
          nodeMap[3] = 4;
          nodeMap[4] = 3;
          nodeMap[5] = 7;
          nodeMap[6] = 6;
          nodeMap[7] = 5;
          nodeMap[8] = 9;
          nodeMap[9] = 8;
          nodeMap[10] = 11;
          nodeMap[11] = 10;
          nodeMap[12] = 14;
          nodeMap[13] = 13;
          nodeMap[14] = 12;
          nodeMap[15] = 16;
          nodeMap[16] = 15;
          nodeMap[17] = 19;
          nodeMap[18] = 18;
          nodeMap[19] = 17;
        }
      }
      else if (dir == 1)
      {
        if (polyOrder == 1)
        {
          nodeMap[0] = 2;
          nodeMap[1] = 3;
          nodeMap[2] = 0;
          nodeMap[3] = 1;
          nodeMap[4] = 6;
          nodeMap[5] = 7;
          nodeMap[6] = 4;
          nodeMap[7] = 5;
        }
        else if (polyOrder == 2)
        {
          nodeMap[0] = 5;
          nodeMap[1] = 6;
          nodeMap[2] = 7;
          nodeMap[3] = 3;
          nodeMap[4] = 4;
          nodeMap[5] = 0;
          nodeMap[6] = 1;
          nodeMap[7] = 2;
          nodeMap[8] = 10;
          nodeMap[9] = 11;
          nodeMap[10] = 8;
          nodeMap[11] = 9;
          nodeMap[12] = 17;
          nodeMap[13] = 18;
          nodeMap[14] = 19;
          nodeMap[15] = 15;
          nodeMap[16] = 16;
          nodeMap[17] = 12;
          nodeMap[18] = 13;
          nodeMap[19] = 14;
        }
      }
    }*/
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
    if (NDIM == 2)
    {
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
        data[2] = fldPtr_t[0];
        data[3] = fldPtr_tr[0];
      }
      else if (polyOrder == 2)
      {
        data[0] = fldPtr[0];
        data[1] = fldPtr[1];
        data[2] = fldPtr_r[0];
        data[3] = fldPtr[2];
        data[4] = fldPtr_r[2];
        data[5] = fldPtr_t[0];
        data[6] = fldPtr_t[1];
        data[7] = fldPtr_tr[0];
      }
      else
      {
        Lucee::Except lce("SerendipityElement::extractFromField: polyOrder ");
        lce << polyOrder << " not supported";
        throw lce;
      }
    }
    else if (NDIM == 3)
    {
      Lucee::ConstFieldPtr<double> fldPtr_000 = fld.createConstPtr();
      Lucee::ConstFieldPtr<double> fldPtr_100 = fld.createConstPtr();
      Lucee::ConstFieldPtr<double> fldPtr_010 = fld.createConstPtr();
      Lucee::ConstFieldPtr<double> fldPtr_110 = fld.createConstPtr();
      Lucee::ConstFieldPtr<double> fldPtr_001 = fld.createConstPtr();
      Lucee::ConstFieldPtr<double> fldPtr_101 = fld.createConstPtr();
      Lucee::ConstFieldPtr<double> fldPtr_011 = fld.createConstPtr();
      Lucee::ConstFieldPtr<double> fldPtr_111 = fld.createConstPtr();
      // attach pointers to proper locations
      fld.setPtr(fldPtr_000, this->currIdx[0], this->currIdx[1], this->currIdx[2]);
      fld.setPtr(fldPtr_100, this->currIdx[0]+1, this->currIdx[1], this->currIdx[2]);
      fld.setPtr(fldPtr_010, this->currIdx[0], this->currIdx[1]+1, this->currIdx[2]);
      fld.setPtr(fldPtr_110, this->currIdx[0]+1, this->currIdx[1]+1, this->currIdx[2]);
      fld.setPtr(fldPtr_001, this->currIdx[0], this->currIdx[1], this->currIdx[2]+1);
      fld.setPtr(fldPtr_101, this->currIdx[0]+1, this->currIdx[1], this->currIdx[2]+1);
      fld.setPtr(fldPtr_011, this->currIdx[0], this->currIdx[1]+1, this->currIdx[2]+1);
      fld.setPtr(fldPtr_111, this->currIdx[0]+1, this->currIdx[1]+1, this->currIdx[2]+1);
      // copy data
      if (polyOrder == 1)
      {
        data[0] = fldPtr_000[0];
        data[1] = fldPtr_100[0];
        data[2] = fldPtr_010[0];
        data[3] = fldPtr_110[0];
        data[4] = fldPtr_001[0];
        data[5] = fldPtr_101[0];
        data[6] = fldPtr_011[0];
        data[7] = fldPtr_111[0];
      }
      else
      {
        Lucee::Except lce("SerendipityElement::extractFromField: polyOrder ");
        lce << polyOrder << " not supported";
        throw lce;
      }
    }
    else
    {
      Lucee::Except lce("SerendipityElement::extractFromField: NDIM ");
      lce << NDIM << " not supported";
      throw lce;
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::copyAllDataFromField(const Lucee::Field<NDIM, double>& fld,
    double *data)
  {
    if (NDIM == 2)
    {
      // region to copy
      Lucee::Region<2, int> rgn =
        this->template getGrid<Lucee::StructuredGridBase<2> >().getLocalRegion();

      Lucee::ConstFieldPtr<double> fldPtr = fld.createConstPtr();

      if (polyOrder == 1)
      {
        unsigned count = 0;
        for (int i=rgn.getLower(0); i<rgn.getUpper(0)+1; ++i)
          for (int j=rgn.getLower(1); j<rgn.getUpper(1)+1; ++j)
          {
            fld.setPtr(fldPtr, i, j);
            data[count++] = fldPtr[0];
          }
      }
      else if (polyOrder == 2)
      {
        std::vector<int> glob, loc;
        for (int i = rgn.getLower(0); i < rgn.getUpper(0)+1; i++)
          for (int j = rgn.getLower(1); j < rgn.getUpper(1)+1; j++)
          {
            fld.setPtr(fldPtr, i, j);
            // determine mapping of exclusively owned nodes
            getGlobalIndices(i, j, glob, loc);
            for (unsigned n=0; n<glob.size(); ++n)
            {
              data[glob[n]] = fldPtr[loc[n]];
            }
          }
      }
      else throw Lucee::Except("SerendipityElement::copyAllDataFromField: polyOrder implemented!");
    }
    else if (NDIM == 3)
    {
      // region to copy
      Lucee::Region<3, int> rgn =
        this->template getGrid<Lucee::StructuredGridBase<3> >().getLocalRegion();

      Lucee::ConstFieldPtr<double> fldPtr = fld.createConstPtr();

      if (polyOrder == 1)
      {
        unsigned count = 0;
        for (int i=rgn.getLower(0); i<rgn.getUpper(0)+1; ++i)
          for (int j=rgn.getLower(1); j<rgn.getUpper(1)+1; ++j)
            for (int k=rgn.getLower(2); k<rgn.getUpper(2)+1; ++k)
            {
              fld.setPtr(fldPtr, i, j, k);
              data[count++] = fldPtr[0];
            }
      }
      else throw Lucee::Except("SerendipityElement::copyAllDataFromField: polyOrder implemented!");
    }
    else throw Lucee::Except("SerendipityElement::copyAllDataFromField: Dimension not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::copyAllDataToField(const double *data, 
    Lucee::Field<NDIM, double>& fld)
  {
    if (NDIM == 2)
    {
      // region to copy
      Lucee::Region<2, int> rgn =
        this->template getGrid<Lucee::StructuredGridBase<2> >().getLocalRegion();

      Lucee::FieldPtr<double> fldPtr = fld.createPtr();

      if (polyOrder == 1)
      {
        unsigned count = 0;
        for (int i=rgn.getLower(0); i<rgn.getUpper(0)+1; ++i)
          for (int j=rgn.getLower(1); j<rgn.getUpper(1)+1; ++j)
          {
            fld.setPtr(fldPtr, i, j);
            fldPtr[0] = data[count++];
          }
      }
      else if (polyOrder == 2)
      {
        std::vector<int> glob, loc;
        for (int i=rgn.getLower(0); i<rgn.getUpper(0)+1; ++i)
          for (int j=rgn.getLower(1); j<rgn.getUpper(1)+1; ++j)
          {
            fld.setPtr(fldPtr, i, j);
            // determine mapping of exclusively owned nodes
            getGlobalIndices(i, j, glob, loc);
            for (unsigned n=0; n<glob.size(); ++n)
            {
              fldPtr[loc[n]] = data[glob[n]];
            }
          }
      }
      else throw Lucee::Except("SerendipityElement::copyAllDataToField: polyOrder not implemented!");
    }
    else if (NDIM == 3)
    {
      // region to copy
      Lucee::Region<3, int> rgn =
        this->template getGrid<Lucee::StructuredGridBase<3> >().getLocalRegion();

      Lucee::FieldPtr<double> fldPtr = fld.createPtr();

      if (polyOrder == 1)
      {
        unsigned count = 0;
        for (int i=rgn.getLower(0); i<rgn.getUpper(0)+1; ++i)
          for (int j=rgn.getLower(1); j<rgn.getUpper(1)+1; ++j)
            for (int k=rgn.getLower(2); k<rgn.getUpper(2)+1; ++k)
            {
              fld.setPtr(fldPtr, i, j, k);
              fldPtr[0] = data[count++];
            }
      }
    }
    else throw Lucee::Except("SerendipityElement::copyAllDataToField: Dimension not implemented!");
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
  SerendipityElement<NDIM>::evalBasis(double xc[NDIM], std::vector<double>& vals,
    std::vector<int>& nodeNum) const
  {
    // Create vector from xc
    Eigen::VectorXd nodeVec(NDIM);
    for (int componentIndex = 0; componentIndex < NDIM; componentIndex++)
      nodeVec(componentIndex) = xc[componentIndex];

    for (int basisIndex = 0; basisIndex < nodeNum.size(); basisIndex++)
      vals[basisIndex] = evalPolynomial(functionVector[nodeNum[basisIndex]], nodeVec);
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::setupMatrices()
  {
    clock_t tmStart, tmEnd;
    tmStart = clock();
    
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

      Lucee::ColMajorSequencer<NDIM-1> surfSeq = ColMajorSequencer<NDIM-1>(surfRegion);
      Lucee::ColMajorIndexer<NDIM-1> surfIdxr = ColMajorIndexer<NDIM-1>(surfRegion);

      
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

    tmEnd = clock();
    TELL (std::cout << "Phase I took " << getElapsedTime(tmStart, tmEnd) << std::endl;);
    tmStart = clock();

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

    tmEnd = clock();
    TELL (std::cout << "Phase II took " << getElapsedTime(tmStart, tmEnd) << std::endl;);
    tmStart = clock();    

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

    clock_t t1, t2, t3, t4; double totalTm = 0.0, totalTmD = 0.0;
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
        t1 = clock();
        functionEvaluations(nodeNumber, basisIndex) = evalPolynomial(functionVector[basisIndex],gaussNodeVec);
        t2 = clock();
        totalTm += getElapsedTime(t1, t2);

        t3 = clock();        
        for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
          functionDEvaluations(nodeNumber, basisIndex, dimIndex) = evalPolynomial(computePolynomialDerivative(functionVector[basisIndex],dimIndex),gaussNodeVec);
        t4 = clock();
        totalTmD += getElapsedTime(t3, t4);
      }
    }

    tmEnd = clock();
    TELL (std::cout << "Phase III took " << getElapsedTime(tmStart, tmEnd) << std::endl;);
    TELL (std::cout << "... with inner loop taking " << totalTm << " " << totalTmD << std::endl;);
    tmStart = clock();        

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

      // Set up face-to-interior mapping matrix
      lowerFaceToInteriorMapMatrices = std::vector<Eigen::MatrixXd>(NDIM);
      
      for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
        lowerFaceToInteriorMapMatrices[dimIndex] = Eigen::MatrixXd::Zero(functionVector.size(), polyOrder + 1);

      // Explicitly assign values to matrix elements (somewhat long)
      #include <LcSerendipityElement2DFaceToInteriorOutput>
    }
    else if (NDIM == 3)
    {
      // Set up diffusion matrices
      iMatDiffusion          = std::vector<Eigen::MatrixXd>(NDIM);
      lowerMatDiffusion      = std::vector<Eigen::MatrixXd>(NDIM);
      upperMatDiffusion      = std::vector<Eigen::MatrixXd>(NDIM);

      for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
      {
        iMatDiffusion[dimIndex]          = Eigen::MatrixXd::Zero(functionVector.size(),functionVector.size());
        lowerMatDiffusion[dimIndex]      = Eigen::MatrixXd::Zero(functionVector.size(),functionVector.size());
        upperMatDiffusion[dimIndex]      = Eigen::MatrixXd::Zero(functionVector.size(),functionVector.size());
      }

      // Explicity assign values to matrix elements (very long)
      #include <LcSerendipityElement3DDiffusionOutput>

      // Set up face-to-interior mapping matrix
      lowerFaceToInteriorMapMatrices = std::vector<Eigen::MatrixXd>(NDIM);
      
      // Assumes that all faces have same number of nodes
      for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
        lowerFaceToInteriorMapMatrices[dimIndex] = Eigen::MatrixXd::Zero(functionVector.size(), getNumSurfLowerNodes(0));

      // Explicitly assign values to matrix elements (somewhat long)
      #include <LcSerendipityElement3DFaceToInteriorOutput>
    }

    computeMass(refMass);
    computeStiffness(functionDEvaluations, refStiffness);
    computePerpStiffness(functionDEvaluations, refPerpStiffness);

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
      refPerpStiffness *= 0.5*dq[dimIndex];

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

    tmEnd = clock();
    TELL (std::cout << "Phase IV took " << getElapsedTime(tmStart, tmEnd) << std::endl;);
    tmStart = clock();        
    //computeTransformationScales();
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::resizeMatrices()
  {
    int numNodes = this->getNumNodes();

    refMass          = Eigen::MatrixXd::Zero(numNodes, numNodes);
    refStiffness     = Eigen::MatrixXd::Zero(numNodes, numNodes);
    refPerpStiffness = Eigen::MatrixXd::Zero(numNodes, numNodes);

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
      {
        nodeMatrix << -1,
                      1;
      }
      else if (degree == 2)
      {
        nodeMatrix << -1,
                       0,
                       1;
      }
      else if (degree == 3)
      {
        nodeMatrix << -1,
                      -1/3.0,
                      1/3.0,
                      1;
      }
      else if (degree == 4)
      {
        nodeMatrix << -1,
                      -0.5,
                      0.5,
                      1;
      }
      else if (degree == 5)
      {
        nodeMatrix << -1,
                      -0.6,
	              -0.2,
	              0.2,
                      0.6,
                      1;
      }
      else if (degree == 6)
      {
        nodeMatrix << -1,
                      -2.0/3.0,
	              -1.0/3.0,
	              0,
	              1.0/3.0,
                      2.0/3.0,
                      1;
      }
      else if (degree == 7)
      {
        nodeMatrix << -1,
                      -5.0/7.0,
	              -3.0/7.0,
	              -1.0/7.0,
	              1.0/7.0,
	              3.0/7.0,
                      5.0/7.0,
                      1;
      }
      else if (degree == 8)
      {
        nodeMatrix << -1,
                      -0.75,
	              -0.50,
	              -0.25,
	              0,
	              0.25,
	              0.50,
                      0.75,
                      1;
      }
      else if (degree == 9)
      {
        nodeMatrix << -1,
                      -7.0/9.0,
	              -5.0/9.0,
	              -3.0/9.0,
	              -1.0/9.0,
	              1.0/9.0,
	              3.0/9.0,
	              5.0/9.0,
                      7.0/9.0,
                      1;
      }
      else if (degree == 10)
      {
        nodeMatrix << -1,
                      -0.8,
	              -0.6,
	              -0.4,
	              -0.2,
	              0,
	              0.2,
	              0.4,
	              0.6,
                      0.8,
                      1;
      }
    }
    else if (dimension == 2)
    {
      if (degree == 1)
      {
        nodeMatrix << -1,-1,
                      1,-1,
                      -1,1,
                      1,1;
      }
      else if (degree == 2)
      {
        nodeMatrix << -1,-1,
                    0,-1,
                    1,-1,
                    -1,0,
                    1,0,
                    -1,1,
                    0,1,
                    1,1;
      }
      else if (degree == 3)
      {
        nodeMatrix << -1,-1,
                    -1/3.0,-1,
                    1/3.0,-1,
                    1,-1,
                    -1,-1/3.0,
                    1,-1/3.0,
                    -1,1/3.0,
                    1,1/3.0,
                    -1,1,
                    -1/3.0,1,
                    1/3.0,1,
                    1,1;
      }
      else if (degree == 4)
      {
        nodeMatrix << -1,-1,
                    -0.5,-1,
                    0,-1,
	            0.5,-1,
                    1,-1,
	            -1,-0.5,
                    1,-0.5,
	            -1,0,
	            0,0,
	            1,0,
                    -1,0.5,
                    1,0.5,
                    -1,1,
                    -0.5,1,
	            0,1,
                    0.5,1,
                    1,1;
      }
      else if (degree == 5)
      {
        nodeMatrix << -1,-1,
                    -0.6,-1,
	            -0.2,-1,
                    0.2,-1,
	            0.6,-1,
                    1,-1,
	            -1,-0.6,
                    1,-0.6,
	            -1,-0.2,
	            -0.2,-0.2,
	            0.2,-0.2,
	            1,-0.2,
	            -1,0.2,
	            0,0.2,
	            1,0.2,
                    -1,0.6,
                    1,0.6,
                    -1,1,
                    -0.6,1,
	            -0.2,1,
	            0.2,1,
                    0.6,1,
                    1,1;
      }
      else if (degree == 6)
      {
        nodeMatrix << -1,-1,
                    -2.0/3.0,-1,
	            -1.0/3.0,-1,
	            0,-1,
                    1.0/3.0,-1,
	            2.0/3.0,-1,
                    1,-1,
	            -1,-2.0/3.0,
                    1,-2.0/3.0,
	            -1,-1.0/3.0,
	            -1.0/3.0,-1.0/3.0,
	            0,-1.0/3.0,
	            1.0/3.0,-1.0/3.0,
	            1,-1.0/3.0,
	            -1,0,
	            -1.0/6.0,0,
	            1.0/6.0,0,
	            1,0,
                    -1,1.0/3.0,
	            0,1.0/3.0,
                    1,1.0/3.0,
	            -1,2.0/3.0,
	            1,2.0/3.0,
                    -1,1,
                    -2.0/3.0,1,
	            -1.0/3.0,1,
	            0,1,
	            1.0/3.0,1,
                    2.0/3.0,1,
                    1,1;
      }
      else if (degree == 7)
      {
        nodeMatrix << -1,-1,
                    -5.0/7.0,-1,
	            -3.0/7.0,-1,
	            -1.0/7.0,-1,
	            1.0/7.0,-1,
                    3.0/7.0,-1,
	            5.0/7.0,-1,
                    1,-1,
	            -1,-5.0/7.0,
                    1,-5.0/7.0,
	            -1,-3.0/7.0,
	            -3.0/7.0,-3.0/7.0,
	            -1.0/7.0,-3.0/7.0,
	            1.0/7.0,-3.0/7.0,
	            3.0/7.0,-3.0/7.0,
	            1,-3.0/7.0,
	            -1,-1.0/7.0,
	            -2.0/7.0,-1.0/7.0,
	            0,-1.0/7.0,
	            2.0/7.0,-1.0/7.0,
	            1,-1.0/7.0,
	            -1,1.0/7.0,
	            -1.0/7.0,1.0/7.0,
	            1.0/7.0,1.0/7.0,
	            1,1.0/7.0,
	            -1,3.0/7.0,
	            0,3.0/7.0,
	            1,3.0/7.0,
	            -1,5.0/7.0,
	            1,5.0/7.0,
                    -1,1,
                    -5.0/7.0,1,
	            -3.0/7.0,1,
	            -1.0/7.0,1,
	            1.0/7.0,1,
	            3.0/7.0,1,
                    5.0/7.0,1,
                    1,1;
      }
      else if (degree == 8)
      {
        nodeMatrix << -1,-1,
                    -0.75,-1,
	            -0.5,-1,
	            -0.25,-1,
                    0,-1,
	            0.25,-1,
	            0.5,-1,
	            0.75,-1,
                    1,-1,
	            -1,-0.75,
	            1,-0.75,
	            -1,-0.5,
	            -0.5,-0.5,
	            -0.25,-0.5,
	            0,-0.5,
	            0.25,-0.5,
	            0.5,-0.5,
                    1,-0.5,
	            -1,-0.25,
	            -0.375,-0.25,
	            -0.125,-0.25,
	            0.125,-0.25,
	            0.375,-0.25,
	            1,-0.25,
	            -1,0,
	            -0.25,0,
	            0,0,
	            0.25,0,
	            1,0,
	            -1,0.25,
	            -0.125,0.25,
	            0.125,0.25,
	            1,0.25,
                    -1,0.5,
	            0,0.5,
                    1,0.5,
	            -1,0.75,
	            1,0.75,
                    -1,1,
	            -0.75,1,
                    -0.5,1,
	            -0.25,1,
	            0,1,
	            0.25,1,
                    0.5,1,
	            0.75,1,
                    1,1;
      }
      else if (degree == 9)
      {
        nodeMatrix << -1,-1,
                    -7.0/9.0,-1,
	            -5.0/9.0,-1,
	            -3.0/9.0,-1,
	            -1.0/9.0,-1,
	            1.0/9.0,-1,
                    3.0/9.0,-1,
	            5.0/9.0,-1,
	            7.0/9.0,-1,
                    1,-1,
	            -1,-7.0/9.0,
	            1,-7.0/9.0,
	            -1,-5.0/9.0,
	            -5.0/9.0,-5.0/9.0,
	            -3.0/9.0,-5.0/9.0,
	            -1.0/9.0,-5.0/9.0,
                    1.0/9.0,-5.0/9.0,
                    3.0/9.0,-5.0/9.0,
                    5.0/9.0,-5.0/9.0,
                    1,-5.0/9.0,
	            -1,-3.0/9.0,
	            -4.0/9.0,-3.0/9.0,
	            -2.0/9.0,-3.0/9.0,
	            0,-3.0/9.0,
	            2.0/9.0,-3.0/9.0,
	            4.0/9.0,-3.0/9.0,
	            1,-3.0/9.0,
	            -1,-1.0/9.0,
	            -3.0/9.0,-1.0/9.0,
	            -1.0/9.0,-1.0/9.0,
	            1.0/9.0,-1.0/9.0,
	            3.0/9.0,-1.0/9.0,
	            1,-1.0/9.0,
	            -1,1.0/9.0,
	            -2.0/9.0,1.0/9.0,
	            0,1.0/9.0,
	            2.0/9.0,1.0/9.0,
	            1,1.0/9.0,
	            -1,3.0/9.0,
	            -1.0/9.0,3.0/9.0,
                    1.0/9.0,3.0/9.0,
	            1,3.0/9.0,
	            -1,5.0/9.0,
	            0,5.0/9.0,
	            1,5.0/9.0,
	            -1,7.0/9.0,
                    1,7.0/9.0,
                    -1,1,
                    -7.0/9.0,1,
	            -5.0/9.0,1,
	            -3.0/9.0,1,
	            -1.0/9.0,1,
	            1.0/9.0,1,
                    3.0/9.0,1,
	            5.0/9.0,1,
	            7.0/9.0,1,
                    1,1;
      }
      else if (degree == 10)
      {
        nodeMatrix << -1,-1,
                    -0.8,-1,
	            -0.6,-1,
	            -0.4,-1,
	            -0.2,-1,
                    0,-1,
	            0.2,-1,
	            0.4,-1,
	            0.6,-1,
	            0.8,-1,
                    1,-1,
	            -1,-0.8,
	            1,-0.8,
	            -1,-0.6,
	            -0.6,-0.6,
	            -0.4,-0.6,
	            -0.2,-0.6,
	            0,-0.6,
	            0.2,-0.6,
	            0.4,-0.6,
	            0.6,-0.6,
                    1,-0.6,
	            -1,-0.4,
	            -0.5,-0.4,
	            -0.3,-0.4,
	            -0.1,-0.4,
	            0.1,-0.4,
	            0.3,-0.4,
	            0.5,-0.4,
	            1,-0.4,
	            -1,-0.2,
	            -0.4,-0.2,
	            -0.2,-0.2,
	            0,-0.2,
                    0.2,-0.2,
                    0.4,-0.2,
	            1,-0.2,
	            -1,0,
	            -0.3,0,
	            -0.1,0,
	            0.1,0,
	            0.3,0,
	            1,0,
	            -1,0.2,
	            -0.2,0.2,
	            0,0.2,
	            0.2,0.2,
	            1,0.2,
                    -1,0.4,
                    -0.1,0.4,
	            0.1,0.4,
                    1,0.4,
	            -1,0.6,
	            0,0.6,
	            1,0.6,
	            -1,0.8,
	            1,0.8,
                    -1,1,
	            -0.8,1,
	            -0.6,1,
	            -0.4,1,
	            -0.2,1,
                    0,1,
	            0.2,1,
	            0.4,1,
	            0.6,1,
	            0.8,1,
                    1,1;
      }
    }
    else if (dimension == 3)
    {
      if (degree == 1)
      {
        nodeMatrix << -1,-1,-1,
                     1,-1,-1,
                     -1,1,-1,
                     1,1,-1,
                     -1,-1,1,
                     1,-1,1,
                     -1,1,1,
                     1,1,1;
      }
      else if (degree == 2)
      {
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
      else if (degree == 3)
      {
        nodeMatrix << -1,-1,-1,
                    -1.0/3.0,-1,-1,
                    1.0/3.0,-1,-1,
                    1,-1,-1,
                    -1,-1.0/3.0,-1,
                    1,-1.0/3.0,-1,
                    -1,1.0/3.0,-1,
                    1,1.0/3.0,-1,
                    -1,1,-1,
                    -1.0/3.0,1,-1,
                    1.0/3.0,1,-1,
                    1,1,-1,
                    -1,-1,-1.0/3.0,
                    1,-1,-1.0/3.0,
                    -1,1,-1.0/3.0,
                    1,1,-1.0/3.0,
                    -1,-1,1.0/3.0,
                    1,-1,1.0/3.0,
                    -1,1,1.0/3.0,
                    1,1,1.0/3.0,
                    -1,-1,1,
                    -1.0/3.0,-1,1,
                    1.0/3.0,-1,1,
                    1,-1,1,
                    -1,-1.0/3.0,1,
                    1,-1.0/3.0,1,
                    -1,1.0/3.0,1,
                    1,1.0/3.0,1,
                    -1,1,1,
                    -1.0/3.0,1,1,
                    1.0/3.0,1,1,
                    1,1,1;
      }
      else if (degree == 4)
      {
        nodeMatrix << -1,-1,-1,
	            -0.5,-1,-1,
                    0,-1,-1,
                    0.5,-1,-1,
                    1,-1,-1,
                    -1,-0.5,-1,
                    1,-0.5,-1,
                    -1,0,-1,
	            0,0,-1,
                    1,0,-1,
                    -1,0.5,-1,
                    1,0.5,-1,
                    -1,1,-1,
	            -0.5,1,-1,
                    0,1,-1,
                    0.5,1,-1,
                    1,1,-1,
                    -1,-1,-0.5,
                    1,-1,-0.5,
                    -1,1,-0.5,
                    1,1,-0.5,
                    -1,-1,0,
	            0,-1,0,
                    1,-1,0,
                    -1,0,0,
                    1,0,0,
                    -1,1,0,
	            0,1,0,
                    1,1,0,
                    -1,-1,0.5,
                    1,-1,0.5,
                    -1,1,0.5,
                    1,1,0.5,
                    -1,-1,1,
	            -0.5,-1,1,
                    0,-1,1,
                    0.5,-1,1,
                    1,-1,1,
                    -1,-0.5,1,
                    1,-0.5,1,
                    -1,0,1,
	            0,0,1,
                    1,0,1,
                    -1,0.5,1,
                    1,0.5,1,
                    -1,1,1,
	            -0.5,1,1,
                    0,1,1,
                    0.5,1,1,
                    1,1,1;
      }
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
      else if (degree == 3)
      {
        nodeMatrix << -1,-1,-1,-1,
                    -1.0/3.0,-1,-1,-1,
                    1.0/3.0,-1,-1,-1,
                    1,-1,-1,-1,
                    -1,-1.0/3.0,-1,-1,
                    1,-1.0/3.0,-1,-1,
                    -1,1.0/3.0,-1,-1,
                    1,1.0/3.0,-1,-1,
                    -1,1,-1,-1,
                    -1.0/3.0,1,-1,-1,
                    1.0/3.0,1,-1,-1,
                    1,1,-1,-1,
                    -1,-1,-1.0/3.0,-1,
                    1,-1,-1.0/3.0,-1,
                    -1,1,-1.0/3.0,-1,
                    1,1,-1.0/3.0,-1,
                    -1,-1,1.0/3.0,-1,
                    1,-1,1.0/3.0,-1,
                    -1,1,1.0/3.0,-1,
                    1,1,1.0/3.0,-1,
                    -1,-1,1,-1,
                    -1.0/3.0,-1,1,-1,
                    1.0/3.0,-1,1,-1,
                    1,-1,1,-1,
                    -1,-1.0/3.0,1,-1,
                    1,-1.0/3.0,1,-1,
                    -1,1.0/3.0,1,-1,
                    1,1.0/3.0,1,-1,
                    -1,1,1,-1,
                    -1.0/3.0,1,1,-1,
                    1.0/3.0,1,1,-1,
                    1,1,1,-1,
                    -1,-1,-1,-1.0/3.0,
                    1,-1,-1,-1.0/3.0,
                    -1,1,-1,-1.0/3.0,
                    1,1,-1,-1.0/3.0,
                    -1,-1,1,-1.0/3.0,
                    1,-1,1,-1.0/3.0,
                    -1,1,1,-1.0/3.0,
                    1,1,1,-1.0/3.0,
                    -1,-1,-1,1.0/3.0,
                    1,-1,-1,1.0/3.0,
                    -1,1,-1,1.0/3.0,
                    1,1,-1,1.0/3.0,
                    -1,-1,1,1.0/3.0,
                    1,-1,1,1.0/3.0,
                    -1,1,1,1.0/3.0,
                    1,1,1,1.0/3.0,
                    -1,-1,-1,1,
                    -1.0/3.0,-1,-1,1,
                    1.0/3.0,-1,-1,1,
                    1,-1,-1,1,
                    -1,-1.0/3.0,-1,1,
                    1,-1.0/3.0,-1,1,
                    -1,1.0/3.0,-1,1,
                    1,1.0/3.0,-1,1,
                    -1,1,-1,1,
                    -1.0/3.0,1,-1,1,
                    1.0/3.0,1,-1,1,
                    1,1,-1,1,
                    -1,-1,-1.0/3.0,1,
                    1,-1,-1.0/3.0,1,
                    -1,1,-1.0/3.0,1,
                    1,1,-1.0/3.0,1,
                    -1,-1,1.0/3.0,1,
                    1,-1,1.0/3.0,1,
                    -1,1,1.0/3.0,1,
                    1,1,1.0/3.0,1,
                    -1,-1,1,1,
                    -1.0/3.0,-1,1,1,
                    1.0/3.0,-1,1,1,
                    1,-1,1,1,
                    -1,-1.0/3.0,1,1,
                    1,-1.0/3.0,1,1,
                    -1,1.0/3.0,1,1,
                    1,1.0/3.0,1,1,
                    -1,1,1,1,
                    -1.0/3.0,1,1,1,
                    1.0/3.0,1,1,1,
                    1,1,1,1;
      }
      else if (degree == 4)
      {
        nodeMatrix << -1,-1,-1,-1,
	            -0.5,-1,-1,-1,
                    0,-1,-1,-1,
                    0.5,-1,-1,-1,
                    1,-1,-1,-1,
                    -1,-0.5,-1,-1,
                    1,-0.5,-1,-1,
                    -1,0,-1,-1,
	            0,0,-1,-1,
                    1,0,-1,-1,
                    -1,0.5,-1,-1,
                    1,0.5,-1,-1,
                    -1,1,-1,-1,
	            -0.5,1,-1,-1,
                    0,1,-1,-1,
                    0.5,1,-1,-1,
                    1,1,-1,-1,
                    -1,-1,-0.5,-1,
                    1,-1,-0.5,-1,
                    -1,1,-0.5,-1,
                    1,1,-0.5,-1,
                    -1,-1,0,-1,
	            0,-1,0,-1,
                    1,-1,0,-1,
                    -1,0,0,-1,
                    1,0,0,-1,
                    -1,1,0,-1,
	            0,1,0,-1,
                    1,1,0,-1,
                    -1,-1,0.5,-1,
                    1,-1,0.5,-1,
                    -1,1,0.5,-1,
                    1,1,0.5,-1,
                    -1,-1,1,-1,
	            -0.5,-1,1,-1,
                    0,-1,1,-1,
                    0.5,-1,1,-1,
                    1,-1,1,-1,
                    -1,-0.5,1,-1,
                    1,-0.5,1,-1,
                    -1,0,1,-1,
	            0,0,1,-1,
                    1,0,1,-1,
                    -1,0.5,1,-1,
                    1,0.5,1,-1,
                    -1,1,1,-1,
	            -0.5,1,1,-1,
                    0,1,1,-1,
                    0.5,1,1,-1,
	            1,1,1,-1,
                    -1,-1,-1,-0.5,
                    1,-1,-1,-0.5,
                    -1,1,-1,-0.5,
                    1,1,-1,-0.5,
                    -1,-1,1,-0.5,
                    1,-1,1,-0.5,
                    -1,1,1,-0.5,
                    1,1,1,-0.5,
                    -1,-1,-1,0,
                    0,-1,-1,0,
                    1,-1,-1,0,
                    -1,0,-1,0,
                    1,0,-1,0,
                    -1,1,-1,0,
                    0,1,-1,0,
                    1,1,-1,0,
                    -1,-1,0,0,
                    1,-1,0,0,
                    -1,1,0,0,
                    1,1,0,0,
                    -1,-1,1,0,
                    0,-1,1,0,
                    1,-1,1,0,
                    -1,0,1,0,
                    1,0,1,0,
                    -1,1,1,0,
                    0,1,1,0,
                    1,1,1,0,
                    -1,-1,-1,0.5,
                    1,-1,-1,0.5,
                    -1,1,-1,0.5,
                    1,1,-1,0.5,
                    -1,-1,1,0.5,
                    1,-1,1,0.5,
                    -1,1,1,0.5,
                    1,1,1,0.5,
                    -1,-1,-1,1,
	            -0.5,-1,-1,1,
                    0,-1,-1,1,
                    0.5,-1,-1,1,
                    1,-1,-1,1,
                    -1,-0.5,-1,1,
                    1,-0.5,-1,1,
                    -1,0,-1,1,
	            0,0,-1,1,
                    1,0,-1,1,
                    -1,0.5,-1,1,
                    1,0.5,-1,1,
                    -1,1,-1,1,
	            -0.5,1,-1,1,
                    0,1,-1,1,
                    0.5,1,-1,1,
                    1,1,-1,1,
                    -1,-1,-0.5,1,
                    1,-1,-0.5,1,
                    -1,1,-0.5,1,
                    1,1,-0.5,1,
                    -1,-1,0,1,
	            0,-1,0,1,
                    1,-1,0,1,
                    -1,0,0,1,
                    1,0,0,1,
                    -1,1,0,1,
	            0,1,0,1,
                    1,1,0,1,
                    -1,-1,0.5,1,
                    1,-1,0.5,1,
                    -1,1,0.5,1,
                    1,1,0.5,1,
                    -1,-1,1,1,
	            -0.5,-1,1,1,
                    0,-1,1,1,
                    0.5,-1,1,1,
                    1,-1,1,1,
                    -1,-0.5,1,1,
                    1,-0.5,1,1,
                    -1,0,1,1,
	            0,0,1,1,
                    1,0,1,1,
                    -1,0.5,1,1,
                    1,0.5,1,1,
                    -1,1,1,1,
	            -0.5,1,1,1,
                    0,1,1,1,
                    0.5,1,1,1,
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
      else if (degree == 2)
      {
        nodeMatrix << -1,-1,-1,-1,-1,
                    0,-1,-1,-1,-1,
                    1,-1,-1,-1,-1,
                    -1,0,-1,-1,-1,
                    1,0,-1,-1,-1,
                    -1,1,-1,-1,-1,
                    0,1,-1,-1,-1,
                    1,1,-1,-1,-1,
                    -1,-1,0,-1,-1,
                    1,-1,0,-1,-1,
                    -1,1,0,-1,-1,
                    1,1,0,-1,-1,
                    -1,-1,1,-1,-1,
                    0,-1,1,-1,-1,
                    1,-1,1,-1,-1,
                    -1,0,1,-1,-1,
                    1,0,1,-1,-1,
                    -1,1,1,-1,-1,
                    0,1,1,-1,-1,
                    1,1,1,-1,-1,
                    -1,-1,-1,0,-1,
                    1,-1,-1,0,-1,
                    -1,1,-1,0,-1,
                    1,1,-1,0,-1,
                    -1,-1,1,0,-1,
                    1,-1,1,0,-1,
                    -1,1,1,0,-1,
                    1,1,1,0,-1,
                    -1,-1,-1,1,-1,
                    0,-1,-1,1,-1,
                    1,-1,-1,1,-1,
                    -1,0,-1,1,-1,
                    1,0,-1,1,-1,
                    -1,1,-1,1,-1,
                    0,1,-1,1,-1,
                    1,1,-1,1,-1,
                    -1,-1,0,1,-1,
                    1,-1,0,1,-1,
                    -1,1,0,1,-1,
                    1,1,0,1,-1,
                    -1,-1,1,1,-1,
                    0,-1,1,1,-1,
                    1,-1,1,1,-1,
                    -1,0,1,1,-1,
                    1,0,1,1,-1,
                    -1,1,1,1,-1,
                    0,1,1,1,-1,
	            1,1,1,1,-1,
	            -1,-1,-1,-1,0,
                    1,-1,-1,-1,0,
                    -1,1,-1,-1,0,
                    1,1,-1,-1,0,
                    -1,-1,1,-1,0,
                    1,-1,1,-1,0,
                    -1,1,1,-1,0,
                    1,1,1,-1,0,
                    -1,-1,-1,1,0,
                    1,-1,-1,1,0,
                    -1,1,-1,1,0,
                    1,1,-1,1,0,
                    -1,-1,1,1,0,
                    1,-1,1,1,0,
                    -1,1,1,1,0,
	            1,1,1,1,0,
                    -1,-1,-1,-1,1,
                    0,-1,-1,-1,1,
                    1,-1,-1,-1,1,
                    -1,0,-1,-1,1,
                    1,0,-1,-1,1,
                    -1,1,-1,-1,1,
                    0,1,-1,-1,1,
                    1,1,-1,-1,1,
                    -1,-1,0,-1,1,
                    1,-1,0,-1,1,
                    -1,1,0,-1,1,
                    1,1,0,-1,1,
                    -1,-1,1,-1,1,
                    0,-1,1,-1,1,
                    1,-1,1,-1,1,
                    -1,0,1,-1,1,
                    1,0,1,-1,1,
                    -1,1,1,-1,1,
                    0,1,1,-1,1,
                    1,1,1,-1,1,
                    -1,-1,-1,0,1,
                    1,-1,-1,0,1,
                    -1,1,-1,0,1,
                    1,1,-1,0,1,
                    -1,-1,1,0,1,
                    1,-1,1,0,1,
                    -1,1,1,0,1,
                    1,1,1,0,1,
                    -1,-1,-1,1,1,
                    0,-1,-1,1,1,
                    1,-1,-1,1,1,
                    -1,0,-1,1,1,
                    1,0,-1,1,1,
                    -1,1,-1,1,1,
                    0,1,-1,1,1,
                    1,1,-1,1,1,
                    -1,-1,0,1,1,
                    1,-1,0,1,1,
                    -1,1,0,1,1,
                    1,1,0,1,1,
                    -1,-1,1,1,1,
                    0,-1,1,1,1,
                    1,-1,1,1,1,
                    -1,0,1,1,1,
                    1,0,1,1,1,
                    -1,1,1,1,1,
                    0,1,1,1,1,
	            1,1,1,1,1;
      }
      else if (degree == 3)
      {
        nodeMatrix << -1,-1,-1,-1,-1,
                    -1.0/3.0,-1,-1,-1,-1,
                    1.0/3.0,-1,-1,-1,-1,
                    1,-1,-1,-1,-1,
                    -1,-1.0/3.0,-1,-1,-1,
                    1,-1.0/3.0,-1,-1,-1,
                    -1,1.0/3.0,-1,-1,-1,
                    1,1.0/3.0,-1,-1,-1,
                    -1,1,-1,-1,-1,
                    -1.0/3.0,1,-1,-1,-1,
                    1.0/3.0,1,-1,-1,-1,
                    1,1,-1,-1,-1,
                    -1,-1,-1.0/3.0,-1,-1,
                    1,-1,-1.0/3.0,-1,-1,
                    -1,1,-1.0/3.0,-1,-1,
                    1,1,-1.0/3.0,-1,-1,
                    -1,-1,1.0/3.0,-1,-1,
                    1,-1,1.0/3.0,-1,-1,
                    -1,1,1.0/3.0,-1,-1,
                    1,1,1.0/3.0,-1,-1,
                    -1,-1,1,-1,-1,
                    -1.0/3.0,-1,1,-1,-1,
                    1.0/3.0,-1,1,-1,-1,
                    1,-1,1,-1,-1,
                    -1,-1.0/3.0,1,-1,-1,
                    1,-1.0/3.0,1,-1,-1,
                    -1,1.0/3.0,1,-1,-1,
                    1,1.0/3.0,1,-1,-1,
                    -1,1,1,-1,-1,
                    -1.0/3.0,1,1,-1,-1,
                    1.0/3.0,1,1,-1,-1,
                    1,1,1,-1,-1,
                    -1,-1,-1,-1.0/3.0,-1,
                    1,-1,-1,-1.0/3.0,-1,
                    -1,1,-1,-1.0/3.0,-1,
                    1,1,-1,-1.0/3.0,-1,
                    -1,-1,1,-1.0/3.0,-1,
                    1,-1,1,-1.0/3.0,-1,
                    -1,1,1,-1.0/3.0,-1,
                    1,1,1,-1.0/3.0,-1,
                    -1,-1,-1,1.0/3.0,-1,
                    1,-1,-1,1.0/3.0,-1,
                    -1,1,-1,1.0/3.0,-1,
                    1,1,-1,1.0/3.0,-1,
                    -1,-1,1,1.0/3.0,-1,
                    1,-1,1,1.0/3.0,-1,
                    -1,1,1,1.0/3.0,-1,
                    1,1,1,1.0/3.0,-1,
                    -1,-1,-1,1,-1,
                    -1.0/3.0,-1,-1,1,-1,
                    1.0/3.0,-1,-1,1,-1,
                    1,-1,-1,1,-1,
                    -1,-1.0/3.0,-1,1,-1,
                    1,-1.0/3.0,-1,1,-1,
                    -1,1.0/3.0,-1,1,-1,
                    1,1.0/3.0,-1,1,-1,
                    -1,1,-1,1,-1,
                    -1.0/3.0,1,-1,1,-1,
                    1.0/3.0,1,-1,1,-1,
                    1,1,-1,1,-1,
                    -1,-1,-1.0/3.0,1,-1,
                    1,-1,-1.0/3.0,1,-1,
                    -1,1,-1.0/3.0,1,-1,
                    1,1,-1.0/3.0,1,-1,
                    -1,-1,1.0/3.0,1,-1,
                    1,-1,1.0/3.0,1,-1,
                    -1,1,1.0/3.0,1,-1,
                    1,1,1.0/3.0,1,-1,
                    -1,-1,1,1,-1,
                    -1.0/3.0,-1,1,1,-1,
                    1.0/3.0,-1,1,1,-1,
                    1,-1,1,1,-1,
                    -1,-1.0/3.0,1,1,-1,
                    1,-1.0/3.0,1,1,-1,
                    -1,1.0/3.0,1,1,-1,
                    1,1.0/3.0,1,1,-1,
                    -1,1,1,1,-1,
                    -1.0/3.0,1,1,1,-1,
                    1.0/3.0,1,1,1,-1,
	            1,1,1,1,-1,
	            -1,-1,-1,-1,-1.0/3.0,
                    1,-1,-1,-1,-1.0/3.0,
                    -1,1,-1,-1,-1.0/3.0,
                    1,1,-1,-1,-1.0/3.0,
                    -1,-1,1,-1,-1.0/3.0,
                    1,-1,1,-1,-1.0/3.0,
                    -1,1,1,-1,-1.0/3.0,
                    1,1,1,-1,-1.0/3.0,
                    -1,-1,-1,1,-1.0/3.0,
                    1,-1,-1,1,-1.0/3.0,
                    -1,1,-1,1,-1.0/3.0,
                    1,1,-1,1,-1.0/3.0,
                    -1,-1,1,1,-1.0/3.0,
                    1,-1,1,1,-1.0/3.0,
                    -1,1,1,1,-1.0/3.0,
	            1,1,1,1,-1.0/3.0,
	            -1,-1,-1,-1,1.0/3.0,
                    1,-1,-1,-1,1.0/3.0,
                    -1,1,-1,-1,1.0/3.0,
                    1,1,-1,-1,1.0/3.0,
                    -1,-1,1,-1,1.0/3.0,
                    1,-1,1,-1,1.0/3.0,
                    -1,1,1,-1,1.0/3.0,
                    1,1,1,-1,1.0/3.0,
                    -1,-1,-1,1,1.0/3.0,
                    1,-1,-1,1,1.0/3.0,
                    -1,1,-1,1,1.0/3.0,
                    1,1,-1,1,1.0/3.0,
                    -1,-1,1,1,1.0/3.0,
                    1,-1,1,1,1.0/3.0,
                    -1,1,1,1,1.0/3.0,
	            1,1,1,1,1.0/3.0,
	            -1,-1,-1,-1,1,
                    -1.0/3.0,-1,-1,-1,1,
                    1.0/3.0,-1,-1,-1,1,
                    1,-1,-1,-1,1,
                    -1,-1.0/3.0,-1,-1,1,
                    1,-1.0/3.0,-1,-1,1,
                    -1,1.0/3.0,-1,-1,1,
                    1,1.0/3.0,-1,-1,1,
                    -1,1,-1,-1,1,
                    -1.0/3.0,1,-1,-1,1,
                    1.0/3.0,1,-1,-1,1,
                    1,1,-1,-1,1,
                    -1,-1,-1.0/3.0,-1,1,
                    1,-1,-1.0/3.0,-1,1,
                    -1,1,-1.0/3.0,-1,1,
                    1,1,-1.0/3.0,-1,1,
                    -1,-1,1.0/3.0,-1,1,
                    1,-1,1.0/3.0,-1,1,
                    -1,1,1.0/3.0,-1,1,
                    1,1,1.0/3.0,-1,1,
                    -1,-1,1,-1,1,
                    -1.0/3.0,-1,1,-1,1,
                    1.0/3.0,-1,1,-1,1,
                    1,-1,1,-1,1,
                    -1,-1.0/3.0,1,-1,1,
                    1,-1.0/3.0,1,-1,1,
                    -1,1.0/3.0,1,-1,1,
                    1,1.0/3.0,1,-1,1,
                    -1,1,1,-1,1,
                    -1.0/3.0,1,1,-1,1,
                    1.0/3.0,1,1,-1,1,
                    1,1,1,-1,1,
                    -1,-1,-1,-1.0/3.0,1,
                    1,-1,-1,-1.0/3.0,1,
                    -1,1,-1,-1.0/3.0,1,
                    1,1,-1,-1.0/3.0,1,
                    -1,-1,1,-1.0/3.0,1,
                    1,-1,1,-1.0/3.0,1,
                    -1,1,1,-1.0/3.0,1,
                    1,1,1,-1.0/3.0,1,
                    -1,-1,-1,1.0/3.0,1,
                    1,-1,-1,1.0/3.0,1,
                    -1,1,-1,1.0/3.0,1,
                    1,1,-1,1.0/3.0,1,
                    -1,-1,1,1.0/3.0,1,
                    1,-1,1,1.0/3.0,1,
                    -1,1,1,1.0/3.0,1,
                    1,1,1,1.0/3.0,1,
                    -1,-1,-1,1,1,
                    -1.0/3.0,-1,-1,1,1,
                    1.0/3.0,-1,-1,1,1,
                    1,-1,-1,1,1,
                    -1,-1.0/3.0,-1,1,1,
                    1,-1.0/3.0,-1,1,1,
                    -1,1.0/3.0,-1,1,1,
                    1,1.0/3.0,-1,1,1,
                    -1,1,-1,1,1,
                    -1.0/3.0,1,-1,1,1,
                    1.0/3.0,1,-1,1,1,
                    1,1,-1,1,1,
                    -1,-1,-1.0/3.0,1,1,
                    1,-1,-1.0/3.0,1,1,
                    -1,1,-1.0/3.0,1,1,
                    1,1,-1.0/3.0,1,1,
                    -1,-1,1.0/3.0,1,1,
                    1,-1,1.0/3.0,1,1,
                    -1,1,1.0/3.0,1,1,
                    1,1,1.0/3.0,1,1,
                    -1,-1,1,1,1,
                    -1.0/3.0,-1,1,1,1,
                    1.0/3.0,-1,1,1,1,
                    1,-1,1,1,1,
                    -1,-1.0/3.0,1,1,1,
                    1,-1.0/3.0,1,1,1,
                    -1,1.0/3.0,1,1,1,
                    1,1.0/3.0,1,1,1,
                    -1,1,1,1,1,
                    -1.0/3.0,1,1,1,1,
                    1.0/3.0,1,1,1,1,
	            1,1,1,1,1;
      }
      else if (degree == 4)
      {
        nodeMatrix << -1,-1,-1,-1,-1,
	            -0.5,-1,-1,-1,-1,
                    0,-1,-1,-1,-1,
                    0.5,-1,-1,-1,-1,
                    1,-1,-1,-1,-1,
                    -1,-0.5,-1,-1,-1,
                    1,-0.5,-1,-1,-1,
                    -1,0,-1,-1,-1,
	            0,0,-1,-1,-1,
                    1,0,-1,-1,-1,
                    -1,0.5,-1,-1,-1,
                    1,0.5,-1,-1,-1,
                    -1,1,-1,-1,-1,
	            -0.5,1,-1,-1,-1,
                    0,1,-1,-1,-1,
                    0.5,1,-1,-1,-1,
                    1,1,-1,-1,-1,
                    -1,-1,-0.5,-1,-1,
                    1,-1,-0.5,-1,-1,
                    -1,1,-0.5,-1,-1,
                    1,1,-0.5,-1,-1,
                    -1,-1,0,-1,-1,
	            0,-1,0,-1,-1,
                    1,-1,0,-1,-1,
                    -1,0,0,-1,-1,
                    1,0,0,-1,-1,
                    -1,1,0,-1,-1,
	            0,1,0,-1,-1,
                    1,1,0,-1,-1,
                    -1,-1,0.5,-1,-1,
                    1,-1,0.5,-1,-1,
                    -1,1,0.5,-1,-1,
                    1,1,0.5,-1,-1,
                    -1,-1,1,-1,-1,
	            -0.5,-1,1,-1,-1,
                    0,-1,1,-1,-1,
                    0.5,-1,1,-1,-1,
                    1,-1,1,-1,-1,
                    -1,-0.5,1,-1,-1,
                    1,-0.5,1,-1,-1,
                    -1,0,1,-1,-1,
	            0,0,1,-1,-1,
                    1,0,1,-1,-1,
                    -1,0.5,1,-1,-1,
                    1,0.5,1,-1,-1,
                    -1,1,1,-1,-1,
	            -0.5,1,1,-1,-1,
                    0,1,1,-1,-1,
                    0.5,1,1,-1,-1,
	            1,1,1,-1,-1,
                    -1,-1,-1,-0.5,-1,
                    1,-1,-1,-0.5,-1,
                    -1,1,-1,-0.5,-1,
                    1,1,-1,-0.5,-1,
                    -1,-1,1,-0.5,-1,
                    1,-1,1,-0.5,-1,
                    -1,1,1,-0.5,-1,
                    1,1,1,-0.5,-1,
                    -1,-1,-1,0,-1,
                    0,-1,-1,0,-1,
                    1,-1,-1,0,-1,
                    -1,0,-1,0,-1,
                    1,0,-1,0,-1,
                    -1,1,-1,0,-1,
                    0,1,-1,0,-1,
                    1,1,-1,0,-1,
                    -1,-1,0,0,-1,
                    1,-1,0,0,-1,
                    -1,1,0,0,-1,
                    1,1,0,0,-1,
                    -1,-1,1,0,-1,
                    0,-1,1,0,-1,
                    1,-1,1,0,-1,
                    -1,0,1,0,-1,
                    1,0,1,0,-1,
                    -1,1,1,0,-1,
                    0,1,1,0,-1,
                    1,1,1,0,-1,
                    -1,-1,-1,0.5,-1,
                    1,-1,-1,0.5,-1,
                    -1,1,-1,0.5,-1,
                    1,1,-1,0.5,-1,
                    -1,-1,1,0.5,-1,
                    1,-1,1,0.5,-1,
                    -1,1,1,0.5,-1,
                    1,1,1,0.5,-1,
                    -1,-1,-1,1,-1,
	            -0.5,-1,-1,1,-1,
                    0,-1,-1,1,-1,
                    0.5,-1,-1,1,-1,
                    1,-1,-1,1,-1,
                    -1,-0.5,-1,1,-1,
                    1,-0.5,-1,1,-1,
                    -1,0,-1,1,-1,
	            0,0,-1,1,-1,
                    1,0,-1,1,-1,
                    -1,0.5,-1,1,-1,
                    1,0.5,-1,1,-1,
                    -1,1,-1,1,-1,
	            -0.5,1,-1,1,-1,
                    0,1,-1,1,-1,
                    0.5,1,-1,1,-1,
                    1,1,-1,1,-1,
                    -1,-1,-0.5,1,-1,
                    1,-1,-0.5,1,-1,
                    -1,1,-0.5,1,-1,
                    1,1,-0.5,1,-1,
                    -1,-1,0,1,-1,
	            0,-1,0,1,-1,
                    1,-1,0,1,-1,
                    -1,0,0,1,-1,
                    1,0,0,1,-1,
                    -1,1,0,1,-1,
	            0,1,0,1,-1,
                    1,1,0,1,-1,
                    -1,-1,0.5,1,-1,
                    1,-1,0.5,1,-1,
                    -1,1,0.5,1,-1,
                    1,1,0.5,1,-1,
                    -1,-1,1,1,-1,
	            -0.5,-1,1,1,-1,
                    0,-1,1,1,-1,
                    0.5,-1,1,1,-1,
                    1,-1,1,1,-1,
                    -1,-0.5,1,1,-1,
                    1,-0.5,1,1,-1,
                    -1,0,1,1,-1,
	            0,0,1,1,-1,
                    1,0,1,1,-1,
                    -1,0.5,1,1,-1,
                    1,0.5,1,1,-1,
                    -1,1,1,1,-1,
	            -0.5,1,1,1,-1,
                    0,1,1,1,-1,
                    0.5,1,1,1,-1,
	            1,1,1,1,-1,
	            -1,-1,-1,-1,-0.5,
                    1,-1,-1,-1,-0.5,
                    -1,1,-1,-1,-0.5,
                    1,1,-1,-1,-0.5,
                    -1,-1,1,-1,-0.5,
                    1,-1,1,-1,-0.5,
                    -1,1,1,-1,-0.5,
                    1,1,1,-1,-0.5,
                    -1,-1,-1,1,-0.5,
                    1,-1,-1,1,-0.5,
                    -1,1,-1,1,-0.5,
                    1,1,-1,1,-0.5,
                    -1,-1,1,1,-0.5,
                    1,-1,1,1,-0.5,
                    -1,1,1,1,-0.5,
	            1,1,1,1,-0.5,
	            -1,-1,-1,-1,0,
	            0,-1,-1,-1,0,
                    1,-1,-1,-1,0,
                    -1,0,-1,-1,0,
                    1,0,-1,-1,0,
                    -1,1,-1,-1,0,
                    0,1,-1,-1,0,
                    1,1,-1,-1,0,
                    -1,-1,0,-1,0,
                    1,-1,0,-1,0,
                    -1,1,0,-1,0,
                    1,1,0,-1,0,
                    -1,-1,1,-1,0,
                    0,-1,1,-1,0,
                    1,-1,1,-1,0,
                    -1,0,1,-1,0,
                    1,0,1,-1,0,
                    -1,1,1,-1,0,
                    0,1,1,-1,0,
                    1,1,1,-1,0,
                    -1,-1,-1,0,0,
                     1,-1,-1,0,0,
                     -1,1,-1,0,0,
                     1,1,-1,0,0,
                     -1,-1,1,0,0,
                     1,-1,1,0,0,
                     -1,1,1,0,0,
	            1,1,1,0,0,
	            -1,-1,-1,1,0,
	            0,-1,-1,1,0,
                    1,-1,-1,1,0,
                    -1,0,-1,1,0,
                    1,0,-1,1,0,
                    -1,1,-1,1,0,
                    0,1,-1,1,0,
                    1,1,-1,1,0,
                    -1,-1,0,1,0,
                    1,-1,0,1,0,
                    -1,1,0,1,0,
                    1,1,0,1,0,
                    -1,-1,1,1,0,
                    0,-1,1,1,0,
                    1,-1,1,1,0,
                    -1,0,1,1,0,
                    1,0,1,1,0,
                    -1,1,1,1,0,
                    0,1,1,1,0,
                    1,1,1,1,0,
	            -1,-1,-1,-1,0.5,
                    1,-1,-1,-1,0.5,
                    -1,1,-1,-1,0.5,
                    1,1,-1,-1,0.5,
                    -1,-1,1,-1,0.5,
                    1,-1,1,-1,0.5,
                    -1,1,1,-1,0.5,
                    1,1,1,-1,0.5,
                    -1,-1,-1,1,0.5,
                    1,-1,-1,1,0.5,
                    -1,1,-1,1,0.5,
                    1,1,-1,1,0.5,
                    -1,-1,1,1,0.5,
                    1,-1,1,1,0.5,
                    -1,1,1,1,0.5,
	            1,1,1,1,0.5,
                    -1,-1,-1,-1,1,
	            -0.5,-1,-1,-1,1,
                    0,-1,-1,-1,1,
                    0.5,-1,-1,-1,1,
                    1,-1,-1,-1,1,
                    -1,-0.5,-1,-1,1,
                    1,-0.5,-1,-1,1,
                    -1,0,-1,-1,1,
	            0,0,-1,-1,1,
                    1,0,-1,-1,1,
                    -1,0.5,-1,-1,1,
                    1,0.5,-1,-1,1,
                    -1,1,-1,-1,1,
	            -0.5,1,-1,-1,1,
                    0,1,-1,-1,1,
                    0.5,1,-1,-1,1,
                    1,1,-1,-1,1,
                    -1,-1,-0.5,-1,1,
                    1,-1,-0.5,-1,1,
                    -1,1,-0.5,-1,1,
                    1,1,-0.5,-1,1,
                    -1,-1,0,-1,1,
	            0,-1,0,-1,1,
                    1,-1,0,-1,1,
                    -1,0,0,-1,1,
                    1,0,0,-1,1,
                    -1,1,0,-1,1,
	            0,1,0,-1,1,
                    1,1,0,-1,1,
                    -1,-1,0.5,-1,1,
                    1,-1,0.5,-1,1,
                    -1,1,0.5,-1,1,
                    1,1,0.5,-1,1,
                    -1,-1,1,-1,1,
	            -0.5,-1,1,-1,1,
                    0,-1,1,-1,1,
                    0.5,-1,1,-1,1,
                    1,-1,1,-1,1,
                    -1,-0.5,1,-1,1,
                    1,-0.5,1,-1,1,
                    -1,0,1,-1,1,
	            0,0,1,-1,1,
                    1,0,1,-1,1,
                    -1,0.5,1,-1,1,
                    1,0.5,1,-1,1,
                    -1,1,1,-1,1,
	            -0.5,1,1,-1,1,
                    0,1,1,-1,1,
                    0.5,1,1,-1,1,
	            1,1,1,-1,1,
                    -1,-1,-1,-0.5,1,
                    1,-1,-1,-0.5,1,
                    -1,1,-1,-0.5,1,
                    1,1,-1,-0.5,1,
                    -1,-1,1,-0.5,1,
                    1,-1,1,-0.5,1,
                    -1,1,1,-0.5,1,
                    1,1,1,-0.5,1,
                    -1,-1,-1,0,1,
                    0,-1,-1,0,1,
                    1,-1,-1,0,1,
                    -1,0,-1,0,1,
                    1,0,-1,0,1,
                    -1,1,-1,0,1,
                    0,1,-1,0,1,
                    1,1,-1,0,1,
                    -1,-1,0,0,1,
                    1,-1,0,0,1,
                    -1,1,0,0,1,
                    1,1,0,0,1,
                    -1,-1,1,0,1,
                    0,-1,1,0,1,
                    1,-1,1,0,1,
                    -1,0,1,0,1,
                    1,0,1,0,1,
                    -1,1,1,0,1,
                    0,1,1,0,1,
                    1,1,1,0,1,
                    -1,-1,-1,0.5,1,
                    1,-1,-1,0.5,1,
                    -1,1,-1,0.5,1,
                    1,1,-1,0.5,1,
                    -1,-1,1,0.5,1,
                    1,-1,1,0.5,1,
                    -1,1,1,0.5,1,
                    1,1,1,0.5,1,
                    -1,-1,-1,1,1,
	            -0.5,-1,-1,1,1,
                    0,-1,-1,1,1,
                    0.5,-1,-1,1,1,
                    1,-1,-1,1,1,
                    -1,-0.5,-1,1,1,
                    1,-0.5,-1,1,1,
                    -1,0,-1,1,1,
	            0,0,-1,1,1,
                    1,0,-1,1,1,
                    -1,0.5,-1,1,1,
                    1,0.5,-1,1,1,
                    -1,1,-1,1,1,
	            -0.5,1,-1,1,1,
                    0,1,-1,1,1,
                    0.5,1,-1,1,1,
                    1,1,-1,1,1,
                    -1,-1,-0.5,1,1,
                    1,-1,-0.5,1,1,
                    -1,1,-0.5,1,1,
                    1,1,-0.5,1,1,
                    -1,-1,0,1,1,
	            0,-1,0,1,1,
                    1,-1,0,1,1,
                    -1,0,0,1,1,
                    1,0,0,1,1,
                    -1,1,0,1,1,
	            0,1,0,1,1,
                    1,1,0,1,1,
                    -1,-1,0.5,1,1,
                    1,-1,0.5,1,1,
                    -1,1,0.5,1,1,
                    1,1,0.5,1,1,
                    -1,-1,1,1,1,
	            -0.5,-1,1,1,1,
                    0,-1,1,1,1,
                    0.5,-1,1,1,1,
                    1,-1,1,1,1,
                    -1,-0.5,1,1,1,
                    1,-0.5,1,1,1,
                    -1,0,1,1,1,
	            0,0,1,1,1,
                    1,0,1,1,1,
                    -1,0.5,1,1,1,
                    1,0.5,1,1,1,
                    -1,1,1,1,1,
	            -0.5,1,1,1,1,
                    0,1,1,1,1,
                    0.5,1,1,1,1,
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
          for (int powIndex = 0; powIndex < idx[dimIndex]; powIndex++)
            monomialTerm *= nodeCoords(dimIndex);
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
            integrationResult += gaussNodeList(gaussIndex,NDIM)*
              functionDerivative(gaussIndex, kIndex, dimIndex)*
              functionDerivative(gaussIndex, mIndex, dimIndex)*4.0/dq2[dimIndex];
        
        resultMatrix(kIndex,mIndex) = integrationResult;
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::computePerpStiffness(const blitz::Array<double, 3>& functionDerivative,
    Eigen::MatrixXd& resultMatrix)
  {
    for (int kIndex = 0; kIndex < resultMatrix.rows(); kIndex++)
    {
      for (int mIndex = 0; mIndex < resultMatrix.cols(); mIndex++)
      {
        // Reset integration result
        double integrationResult = 0.0;
        
        // Loop over volume gauss points to evaluate integral using gaussian quadrature
        for (int gaussIndex = 0; gaussIndex < gaussNodeList.rows(); gaussIndex++)
          // For now, just assume that perp directions are all but the last direction
          for (int dimIndex = 0; dimIndex < NDIM-1; dimIndex++)
            integrationResult += gaussNodeList(gaussIndex,NDIM)*
              functionDerivative(gaussIndex, kIndex, dimIndex)*
              functionDerivative(gaussIndex, mIndex, dimIndex)*4.0/dq2[dimIndex];
        
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

  template <unsigned NDIM>
  bool
  SerendipityElement<NDIM>::isReflectionNode(int srcIndex, int tarIndex, int reflectDim) const
  {
    // Check to see if all coordinates match, with the reflectDim coordinate reversed in sign
    for (int d = 0; d < NDIM; d++)
    {
      if (reflectDim == d)
      {
        if (std::fabs( -nodeList(srcIndex,d) - nodeList(tarIndex,d)) > 1e-4)
          return false;
      }
      else if (std::fabs( nodeList(srcIndex,d) - nodeList(tarIndex,d)) > 1e-4)
        return false;
    }
    return true;
  }
  
  // instantiations
  template class SerendipityElement<1>;
  template class SerendipityElement<2>;
  template class SerendipityElement<3>;
  template class SerendipityElement<4>;
  template class SerendipityElement<5>;
}
