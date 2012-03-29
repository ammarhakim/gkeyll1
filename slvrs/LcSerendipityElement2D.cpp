/**
 * @file	LcSerendipityElement2D.cpp
 *
 * @brief       Reference finite element with serendipity basis
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSerendipityElement2D.h>
#include <LcStructuredGridBase.h>

namespace Lucee
{
  const char *SerendipityElement2D::id = "Serendipity";

/**
 * This function maps the index (ix,iy) into the global index space in
 * the row having 3 nodes (bottom).
 */
  static
  unsigned F_func(unsigned nx, unsigned ny, int ix, int iy)
  {
    return (2*nx+1)*iy + (nx+1)*iy + 2*ix;
  }

/**
 * This function maps the index (ix,iy) into the global index space in
 * the row having 2 nodes (middle).
 */
  static
  unsigned G_func(unsigned nx, unsigned ny, int ix, int iy)
  {
    return (2*nx+1)*iy + (nx+1)*iy + (2*nx+1) + ix;
  }

  SerendipityElement2D::SerendipityElement2D()
    : Lucee::NodalFiniteElementIfc<2>(4), polyOrder(1), 
      refNjNk(4,4), refDNjDNk(4,4), refDNjNk_0(4,4), refDNjNk_1(4,4), 
      refFaceNjNk_xl(4,2), refFaceNjNk_xu(4,2), refFaceNjNk_yl(4,2), refFaceNjNk_yu(4,2),
      idxr(
        &Lucee::FixedVector<2, unsigned>((unsigned)0)[0], &Lucee::FixedVector<2, int>(1)[0])
  {
// notice funcky initialization of indexer: this is just so code
// compiles as indexers do not have default ctors. It gets reset in
// the readInput() method anyway.
  }

  void
  SerendipityElement2D::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::NodalFiniteElementIfc<2>::readInput(tbl);

// read in polynomial order
    polyOrder = tbl.getNumber("polyOrder");

    if (polyOrder == 1)
      this->setNumNodes(4);
    else if (polyOrder == 2)
      this->setNumNodes(8);
    else
    {
      Lucee::Except lce("SerendipityElement2D: Order must be 1 or 2.");
      lce << " Provided " << polyOrder << " instead";
      throw lce;
    }

// initialize matrices
    if (polyOrder == 1)
      setupPoly1();
    else if (polyOrder == 2)
      setupPoly2();

// determine number of global degrees of freedom
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();
    Lucee::Region<2, int> grgn = grid.getGlobalRegion();
    unsigned nx = grgn.getShape(0), ny = grgn.getShape(1);

    numX = nx; numY = ny;

    if (polyOrder == 1)
      numGlobal = (nx+1)*(ny+1);
    else if (polyOrder ==2)
// there are 3 owned nodes per cell, 2*nx and 2*ny edge nodes along
// the top and right edges and the 1 accounts for the node on the
// top-right corner.
      numGlobal = 3*nx*ny + 2*nx + 2*ny + 1;

// create indexer for local -> global mapping: we are using a
// row-major indexing scheme to be consistent with Petsc default
// matrix layout. This is not strictly needed and even column major
// will work here.
    int shape[2];
    if (polyOrder == 1)
    {
// for 4 nodes per element we can use a simple indexer to map local to
// global. For 8 node element the F_func() and G_func() methods need
// to be used instead.
      shape[0] = nx+1; shape[1] = ny+1;
      Lucee::Region<2, int> lgRgn(shape);
      idxr = RowMajorIndexer<2>(lgRgn);
    }
  }

  void
  SerendipityElement2D::getExclusiveNodeIndices(std::vector<unsigned>& ndIds)
  {
    ndIds.clear();
    if (polyOrder == 1)
    {
      ndIds.resize(1);
      ndIds[0] = 1;
    }
    else if (polyOrder == 2)
    {
      ndIds.resize(3);
      ndIds[0] = 1;
      ndIds[1] = 5;
      ndIds[2] = 8;
    }
  }

  unsigned
  SerendipityElement2D::getNumSurfLowerNodes(unsigned dir) const
  {
    if (polyOrder == 1)
      return 2;
    else if (polyOrder == 2)
      return 3;
  }

  unsigned
  SerendipityElement2D::getNumSurfUpperNodes(unsigned dir) const
  {
    if (polyOrder == 1)
      return 2;
    else if (polyOrder == 2)
      return 3;
  }

  unsigned
  SerendipityElement2D::getNumGlobalNodes() const
  {
    return numGlobal;
  }

  void
  SerendipityElement2D::getLocalToGlobal(std::vector<int>& lgMap) const
  {
    int ix = this->currIdx[0], iy = this->currIdx[1];
    if (polyOrder == 1)
    {
      lgMap[0] = idxr.getIndex(ix, iy); // 1
      lgMap[1] = idxr.getIndex(ix+1, iy); // 2
      lgMap[2] = idxr.getIndex(ix+1, iy+1); // 3
      lgMap[3] = idxr.getIndex(ix, iy+1); // 4
    }
    else if (polyOrder == 2)
    {
      lgMap[0] = F_func(numX, numY, ix, iy); // 1
      lgMap[1] = F_func(numX, numY, ix, iy) + 2; // 2
      lgMap[2] = F_func(numX, numY, ix, iy+1) + 2; // 3
      lgMap[3] = F_func(numX, numY, ix, iy+1); // 4

      lgMap[4] = F_func(numX, numY, ix, iy) + 1; // 5
      lgMap[5] = G_func(numX, numY, ix, iy) + 1; // 6
      lgMap[6] = F_func(numX, numY, ix, iy+1) + 1; // 7
      lgMap[7] = G_func(numX, numY, ix, iy); // 8
    }
  }

  void
  SerendipityElement2D::getSurfLowerLocalToGlobal(unsigned dir,
    std::vector<int>& lgMap) const
  {
    int ix = this->currIdx[0], iy = this->currIdx[1];
    if (polyOrder == 1)
    {
      if (dir == 0)
      {
        lgMap[0] = idxr.getIndex(ix, iy); // 1
        lgMap[1] = idxr.getIndex(ix, iy+1); // 4
      }
      else if (dir == 1)
      {
        lgMap[0] = idxr.getIndex(ix, iy); // 1
        lgMap[1] = idxr.getIndex(ix+1, iy); // 2
      }
    }
    else if (polyOrder == 2)
    {
      if (dir == 0)
      {
        lgMap[0] = F_func(numX, numY, ix, iy); // 1
        lgMap[1] = G_func(numX, numY, ix, iy); // 8
        lgMap[2] = F_func(numX, numY, ix, iy+1); // 4
      }
      else if (dir == 1)
      {
        lgMap[0] = F_func(numX, numY, ix, iy); // 1
        lgMap[1] = F_func(numX, numY, ix, iy) + 1; // 5
        lgMap[2] = F_func(numX, numY, ix, iy) + 2; // 2
      }
    }
  }

  void
  SerendipityElement2D::getSurfUpperLocalToGlobal(unsigned dir,
    std::vector<int>& lgMap) const
  {
    int ix = this->currIdx[0], iy = this->currIdx[1];
    if (polyOrder == 1)
    {
      if (dir == 0)
      {
        lgMap[0] = idxr.getIndex(ix+1, iy); // 2
        lgMap[1] = idxr.getIndex(ix+1, iy+1); // 3
      }
      else if (dir == 1)
      {
        lgMap[0] = idxr.getIndex(ix, iy+1); // 4
        lgMap[1] = idxr.getIndex(ix+1, iy+1); // 3
      }
    }
    else if (polyOrder == 2)
    {
      if (dir == 0)
      {
        lgMap[0] = F_func(numX, numY, ix, iy) + 2; // 2
        lgMap[1] = G_func(numX, numY, ix, iy) + 1; // 6
        lgMap[2] = F_func(numX, numY, ix, iy+1) + 2; // 3
      }
      else if (dir == 1)
      {
        lgMap[0] = F_func(numX, numY, ix, iy+1); // 4
        lgMap[1] = F_func(numX, numY, ix, iy+1) + 1; // 7
        lgMap[2] = F_func(numX, numY, ix, iy+1) + 2; // 3
      }
    }
  }

  void
  SerendipityElement2D::getSurfLowerNodeNums(unsigned dir,
    std::vector<int>& nodeNum) const
  {
    if (polyOrder == 1)
    {
      if (dir == 0)
      {
        nodeNum[0] = 1;
        nodeNum[1] = 4;
      }
      else if (dir == 1)
      {
        nodeNum[0] = 1;
        nodeNum[1] = 2;
      }
    }
    else if (polyOrder == 2)
    {
      if (dir == 0)
      {
        nodeNum[0] = 1;
        nodeNum[1] = 8;
        nodeNum[2] = 4;
      }
      else if (dir == 1)
      {
        nodeNum[0] = 1;
        nodeNum[1] = 5;
        nodeNum[2] = 2;
      }
    }
  }

  void
  SerendipityElement2D::getSurfUpperNodeNums(unsigned dir,
    std::vector<int>& nodeNum) const
  {
    if (polyOrder == 1)
    {
      if (dir == 0)
      {
        nodeNum[0] = 2;
        nodeNum[1] = 3;
      }
      else if (dir == 1)
      {
        nodeNum[0] = 4;
        nodeNum[1] = 3;
      }
    }
    else if (polyOrder == 2)
    {
      if (dir == 0)
      {
        nodeNum[0] = 2;
        nodeNum[1] = 6;
        nodeNum[2] = 3;
      }
      else if (dir == 1)
      {
        nodeNum[0] = 4;
        nodeNum[1] = 7;
        nodeNum[2] = 3;
      }
    }
  }

  void
  SerendipityElement2D::getNodalCoordinates(Lucee::Matrix<double>& nodeCoords)
  {
// get grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();
// set index and get centroid coordinate
    grid.setIndex(this->currIdx);
    double xc[3], dx2, dy2;
    grid.getCentroid(xc);
    dx2 = 0.5*grid.getDx(0);
    dy2 = 0.5*grid.getDx(1);

// first 4 nodes are same in polyOrder 1 and 2
    nodeCoords(0,0) = xc[0]-dx2;
    nodeCoords(0,1) = xc[1]-dy2;
    nodeCoords(0,2) = 0.0;

    nodeCoords(1,0) = xc[0]+dx2;
    nodeCoords(1,1) = xc[1]-dy2;
    nodeCoords(1,2) = 0.0;
    
    nodeCoords(2,0) = xc[0]+dx2;
    nodeCoords(2,1) = xc[1]+dy2;
    nodeCoords(2,2) = 0.0;

    nodeCoords(3,0) = xc[0]-dx2;
    nodeCoords(3,1) = xc[1]+dy2;
    nodeCoords(3,2) = 0.0;

    if (polyOrder == 2)
    {
// edge nodes coordinates
      nodeCoords(4,0) = xc[0];
      nodeCoords(4,1) = xc[1]-dy2;
      nodeCoords(4,2) = 0.0;

      nodeCoords(5,0) = xc[0]+dx2;
      nodeCoords(5,1) = xc[1];
      nodeCoords(5,2) = 0.0;

      nodeCoords(6,0) = xc[0];
      nodeCoords(6,1) = xc[1]+dy2;
      nodeCoords(6,2) = 0.0;

      nodeCoords(7,0) = xc[0]-dx2;
      nodeCoords(7,1) = xc[1];
      nodeCoords(7,2) = 0.0;
    }
  }

  void
  SerendipityElement2D::getWeights(std::vector<double>& w)
  {
    unsigned nn = this->getNumNodes();
    for (unsigned k=0; k<nn; ++k) w[k] = weights[k];
  }

  void
  SerendipityElement2D::getMassMatrix(Lucee::Matrix<double>& NjNk) const
  {
    NjNk.copy(refNjNk);
  }

  void
  SerendipityElement2D::getLowerFaceMassMatrix(unsigned dir,
    Lucee::Matrix<double>& NjNk) const
  {
    if (dir == 0)
      NjNk.copy(refFaceNjNk_xl);
    else if (dir == 1)
      NjNk.copy(refFaceNjNk_yl);
  }

  void
  SerendipityElement2D::getUpperFaceMassMatrix(unsigned dir,
    Lucee::Matrix<double>& NjNk) const
  {
    if (dir == 0)
      NjNk.copy(refFaceNjNk_xu);
    else if (dir == 1)
      NjNk.copy(refFaceNjNk_yu);
  }

  void
  SerendipityElement2D::getStiffnessMatrix(Lucee::Matrix<double>& DNjDNk) const
  {
    DNjDNk.copy(refDNjDNk);
  }

  void
  SerendipityElement2D::getGradStiffnessMatrix(unsigned dir, 
    Lucee::Matrix<double>& DNjNk) const
  {
    if (dir == 0)
      DNjNk.copy(refDNjNk_0);
    else if (dir == 1)
      DNjNk.copy(refDNjNk_1);
    else 
      throw Lucee::Except(
        "SerendipityElement2D::getGradStiffnessMatrix: Can't use this basis in 3D!");
  }

  void
  SerendipityElement2D::extractFromField(const Lucee::Field<2, double>& fld,
    std::vector<double>& data)
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
      data[2] = fldPtr_tr[0];
      data[3] = fldPtr_t[0];
    }
    else if (polyOrder == 2)
    {
      data[0] = fldPtr[0];
      data[1] = fldPtr_r[0];
      data[2] = fldPtr_tr[0];
      data[3] = fldPtr_t[0];

      data[4] = fldPtr[1];
      data[5] = fldPtr_r[2];
      data[6] = fldPtr_t[1];
      data[7] = fldPtr[2];
    }
  }

  void
  SerendipityElement2D::copyAllDataFromField(const Lucee::Field<2, double>& fld, double *data)
  {
// region to copy
    Lucee::Region<2, int> rgn =
      this->getGrid<Lucee::StructuredGridBase<2> >().getLocalRegion();

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
      for (int i=rgn.getLower(0); i<rgn.getUpper(0)+1; ++i)
        for (int j=rgn.getLower(1); j<rgn.getUpper(1)+1; ++j)
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
  }

  void
  SerendipityElement2D::copyAllDataToField(const double *data, Lucee::Field<2, double>& fld)
  {
// region to copy
    Lucee::Region<2, int> rgn =
      this->getGrid<Lucee::StructuredGridBase<2> >().getLocalRegion();

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
  }

  void
  SerendipityElement2D::setupPoly1()
  {
    unsigned shape[2] = {4,4};
    int start[2] = {1,1};

    unsigned faceShape[2] = {4,2};

// get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();

// get grid spacing (this is assumed to be uniform for now)
    double dx = grid.getDx(0), dy = grid.getDx(1);
    double dx2 = dx*dx, dy2 = dy*dy;

// compute weights
    weights.resize(4);
    for (unsigned i=0; i<4; ++i)
      weights[i] = 0.5*dx*0.5*dy*1.0;

// mass matrix (automatically generated. See scripts/serendipity-2D.mac)
    refNjNk = Lucee::Matrix<double>(shape, start);
    refNjNk(1,1) = 4.0/9.0;
    refNjNk(1,2) = 2.0/9.0;
    refNjNk(1,3) = 1.0/9.0;
    refNjNk(1,4) = 2.0/9.0;
    refNjNk(2,1) = 2.0/9.0;
    refNjNk(2,2) = 4.0/9.0;
    refNjNk(2,3) = 2.0/9.0;
    refNjNk(2,4) = 1.0/9.0;
    refNjNk(3,1) = 1.0/9.0;
    refNjNk(3,2) = 2.0/9.0;
    refNjNk(3,3) = 4.0/9.0;
    refNjNk(3,4) = 2.0/9.0;
    refNjNk(4,1) = 2.0/9.0;
    refNjNk(4,2) = 1.0/9.0;
    refNjNk(4,3) = 2.0/9.0;
    refNjNk(4,4) = 4.0/9.0;

// scale to bring this into physical space
    refNjNk *= 0.5*dx*0.5*dy;

// face mass matrix (x-lower)
    refFaceNjNk_xl = Lucee::Matrix<double>(faceShape, start);
    refFaceNjNk_xl(1,1) = 2.0/3.0;
    refFaceNjNk_xl(1,2) = 1.0/3.0;
    refFaceNjNk_xl(2,1) = 0;
    refFaceNjNk_xl(2,2) = 0;
    refFaceNjNk_xl(3,1) = 0;
    refFaceNjNk_xl(3,2) = 0;
    refFaceNjNk_xl(4,1) = 1.0/3.0;
    refFaceNjNk_xl(4,2) = 2.0/3.0;

// scale to bring this into physical space
    refFaceNjNk_xl *= 0.5*dy;

// face mass matrix (x-upper)
    refFaceNjNk_xu = Lucee::Matrix<double>(faceShape, start);
    refFaceNjNk_xu(1,1) = 0;
    refFaceNjNk_xu(1,2) = 0;
    refFaceNjNk_xu(2,1) = 2.0/3.0;
    refFaceNjNk_xu(2,2) = 1.0/3.0;
    refFaceNjNk_xu(3,1) = 1.0/3.0;
    refFaceNjNk_xu(3,2) = 2.0/3.0;
    refFaceNjNk_xu(4,1) = 0;
    refFaceNjNk_xu(4,2) = 0;

// scale to bring this into physical space
    refFaceNjNk_xu *= 0.5*dy;

// face mass matrix (y-lower)
    refFaceNjNk_yl = Lucee::Matrix<double>(faceShape, start);
    refFaceNjNk_yl(1,1) = 2.0/3.0;
    refFaceNjNk_yl(1,2) = 1.0/3.0;
    refFaceNjNk_yl(2,1) = 1.0/3.0;
    refFaceNjNk_yl(2,2) = 2.0/3.0;
    refFaceNjNk_yl(3,1) = 0;
    refFaceNjNk_yl(3,2) = 0;
    refFaceNjNk_yl(4,1) = 0;
    refFaceNjNk_yl(4,2) = 0;

// scale to bring this into physical space
    refFaceNjNk_yl *= 0.5*dx;

// face mass matrix (y-upper)
    refFaceNjNk_yu = Lucee::Matrix<double>(faceShape, start);
    refFaceNjNk_yu(1,1) = 0;
    refFaceNjNk_yu(1,2) = 0;
    refFaceNjNk_yu(2,1) = 0;
    refFaceNjNk_yu(2,2) = 0;
    refFaceNjNk_yu(3,1) = 1.0/3.0;
    refFaceNjNk_yu(3,2) = 2.0/3.0;
    refFaceNjNk_yu(4,1) = 2.0/3.0;
    refFaceNjNk_yu(4,2) = 1.0/3.0;

// scale to bring this into physical space
    refFaceNjNk_yu *= 0.5*dx;

// stiffness matrix (automatically generated. See scripts/serendipity-2D.mac)
    refDNjDNk = Lucee::Matrix<double>(shape, start);
    refDNjDNk(1,1) = (7*dy2+4*dx2)/(dx2*dy2)/6.0+(dy2+4*dx2)/(dx2*dy2)/6.0;
    refDNjDNk(1,2) = -(8*dy2-4*dx2)/(dx2*dy2)/6.0;
    refDNjDNk(1,3) = (-4*dy2-4*dx2)/(dx2*dy2)/6.0;
    refDNjDNk(1,4) = 2.0*(dy2-2*dx2)/(3.0*dx2*dy2);
    refDNjDNk(2,1) = -(8*dy2-4*dx2)/(dx2*dy2)/6.0;
    refDNjDNk(2,2) = (7*dy2+4*dx2)/(dx2*dy2)/6.0+(dy2+4*dx2)/(dx2*dy2)/6.0;
    refDNjDNk(2,3) = 2.0*(dy2-2*dx2)/(3.0*dx2*dy2);
    refDNjDNk(2,4) = (-4*dy2-4*dx2)/(dx2*dy2)/6.0;
    refDNjDNk(3,1) = (-4*dy2-4*dx2)/(dx2*dy2)/6.0;
    refDNjDNk(3,2) = 2.0*(dy2-2*dx2)/(3.0*dx2*dy2);
    refDNjDNk(3,3) = (7*dy2+4*dx2)/(dx2*dy2)/6.0+(dy2+4*dx2)/(dx2*dy2)/6.0;
    refDNjDNk(3,4) = -(8*dy2-4*dx2)/(dx2*dy2)/6.0;
    refDNjDNk(4,1) = 2.0*(dy2-2*dx2)/(3.0*dx2*dy2);
    refDNjDNk(4,2) = (-4*dy2-4*dx2)/(dx2*dy2)/6.0;
    refDNjDNk(4,3) = -(8*dy2-4*dx2)/(dx2*dy2)/6.0;
    refDNjDNk(4,4) = (7*dy2+4*dx2)/(dx2*dy2)/6.0+(dy2+4*dx2)/(dx2*dy2)/6.0;

// scale to bring this into physical space
    refDNjDNk *= 0.5*dx*0.5*dy;

// grad-stiffness matrix (automatically generated. See scripts/serendipity-2D.mac)
    refDNjNk_0 = Lucee::Matrix<double>(shape, start);
    refDNjNk_0(1,1) = (-2.0)/(3.0*dx);
    refDNjNk_0(1,2) = (-2.0)/(3.0*dx);
    refDNjNk_0(1,3) = -1/dx/3.0;
    refDNjNk_0(1,4) = -1/dx/3.0;
    refDNjNk_0(2,1) = 2.0/(3.0*dx);
    refDNjNk_0(2,2) = 2.0/(3.0*dx);
    refDNjNk_0(2,3) = 1/dx/3.0;
    refDNjNk_0(2,4) = 1/dx/3.0;
    refDNjNk_0(3,1) = 1/dx/3.0;
    refDNjNk_0(3,2) = 1/dx/3.0;
    refDNjNk_0(3,3) = 2.0/(3.0*dx);
    refDNjNk_0(3,4) = 2.0/(3.0*dx);
    refDNjNk_0(4,1) = -1/dx/3.0;
    refDNjNk_0(4,2) = -1/dx/3.0;
    refDNjNk_0(4,3) = (-2.0)/(3.0*dx);
    refDNjNk_0(4,4) = (-2.0)/(3.0*dx);

// scale to bring this into physical space
    refDNjNk_0 *= 0.5*dx*0.5*dy;

// grad-stiffness matrix (automatically generated. See scripts/serendipity-2D.mac)
    refDNjNk_1 = Lucee::Matrix<double>(shape, start);
    refDNjNk_1(1,1) = (-2.0)/(3.0*dy);
    refDNjNk_1(1,2) = -1/dy/3.0;
    refDNjNk_1(1,3) = -1/dy/3.0;
    refDNjNk_1(1,4) = (-2.0)/(3.0*dy);
    refDNjNk_1(2,1) = -1/dy/3.0;
    refDNjNk_1(2,2) = (-2.0)/(3.0*dy);
    refDNjNk_1(2,3) = (-2.0)/(3.0*dy);
    refDNjNk_1(2,4) = -1/dy/3.0;
    refDNjNk_1(3,1) = 1/dy/3.0;
    refDNjNk_1(3,2) = 2.0/(3.0*dy);
    refDNjNk_1(3,3) = 2.0/(3.0*dy);
    refDNjNk_1(3,4) = 1/dy/3.0;
    refDNjNk_1(4,1) = 2.0/(3.0*dy);
    refDNjNk_1(4,2) = 1/dy/3.0;
    refDNjNk_1(4,3) = 1/dy/3.0;
    refDNjNk_1(4,4) = 2.0/(3.0*dy);

// scale to bring this into physical space
    refDNjNk_1 *= 0.5*dx*0.5*dy;
  }

  void
  SerendipityElement2D::setupPoly2()
  {
    unsigned shape[2] = {8,8};
    int start[2] = {1,1};

    unsigned faceShape[2] = {8,3};

// get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();

// get grid spacing (this is assumed to be uniform for now)
    double dx = grid.getDx(0), dy = grid.getDx(1);
    double dx2 = dx*dx, dy2 = dy*dy;

// compute weights
    weights.resize(8);
    weights[0] = 0.5*dx*0.5*dy*(-1.0)/3.0;
    weights[1] = 0.5*dx*0.5*dy*(-1.0)/3.0;
    weights[2] = 0.5*dx*0.5*dy*(-1.0)/3.0;
    weights[3] = 0.5*dx*0.5*dy*(-1.0)/3.0;
    weights[4] = 0.5*dx*0.5*dy*4.0/3.0;
    weights[5] = 0.5*dx*0.5*dy*4.0/3.0;
    weights[6] = 0.5*dx*0.5*dy*4.0/3.0;
    weights[7] = 0.5*dx*0.5*dy*4.0/3.0;

// mass matrix (automatically generated. See scripts/serendipity-2D.mac)
    refNjNk = Lucee::Matrix<double>(shape, start);
    refNjNk(1,1) = 2.0/15.0;
    refNjNk(1,2) = 2.0/45.0;
    refNjNk(1,3) = 1.0/15.0;
    refNjNk(1,4) = 2.0/45.0;
    refNjNk(1,5) = (-2.0)/15.0;
    refNjNk(1,6) = (-8.0)/45.0;
    refNjNk(1,7) = (-8.0)/45.0;
    refNjNk(1,8) = (-2.0)/15.0;
    refNjNk(2,1) = 2.0/45.0;
    refNjNk(2,2) = 2.0/15.0;
    refNjNk(2,3) = 2.0/45.0;
    refNjNk(2,4) = 1.0/15.0;
    refNjNk(2,5) = (-2.0)/15.0;
    refNjNk(2,6) = (-2.0)/15.0;
    refNjNk(2,7) = (-8.0)/45.0;
    refNjNk(2,8) = (-8.0)/45.0;
    refNjNk(3,1) = 1.0/15.0;
    refNjNk(3,2) = 2.0/45.0;
    refNjNk(3,3) = 2.0/15.0;
    refNjNk(3,4) = 2.0/45.0;
    refNjNk(3,5) = (-8.0)/45.0;
    refNjNk(3,6) = (-2.0)/15.0;
    refNjNk(3,7) = (-2.0)/15.0;
    refNjNk(3,8) = (-8.0)/45.0;
    refNjNk(4,1) = 2.0/45.0;
    refNjNk(4,2) = 1.0/15.0;
    refNjNk(4,3) = 2.0/45.0;
    refNjNk(4,4) = 2.0/15.0;
    refNjNk(4,5) = (-8.0)/45.0;
    refNjNk(4,6) = (-8.0)/45.0;
    refNjNk(4,7) = (-2.0)/15.0;
    refNjNk(4,8) = (-2.0)/15.0;
    refNjNk(5,1) = (-2.0)/15.0;
    refNjNk(5,2) = (-2.0)/15.0;
    refNjNk(5,3) = (-8.0)/45.0;
    refNjNk(5,4) = (-8.0)/45.0;
    refNjNk(5,5) = 32.0/45.0;
    refNjNk(5,6) = 4.0/9.0;
    refNjNk(5,7) = 16.0/45.0;
    refNjNk(5,8) = 4.0/9.0;
    refNjNk(6,1) = (-8.0)/45.0;
    refNjNk(6,2) = (-2.0)/15.0;
    refNjNk(6,3) = (-2.0)/15.0;
    refNjNk(6,4) = (-8.0)/45.0;
    refNjNk(6,5) = 4.0/9.0;
    refNjNk(6,6) = 32.0/45.0;
    refNjNk(6,7) = 4.0/9.0;
    refNjNk(6,8) = 16.0/45.0;
    refNjNk(7,1) = (-8.0)/45.0;
    refNjNk(7,2) = (-8.0)/45.0;
    refNjNk(7,3) = (-2.0)/15.0;
    refNjNk(7,4) = (-2.0)/15.0;
    refNjNk(7,5) = 16.0/45.0;
    refNjNk(7,6) = 4.0/9.0;
    refNjNk(7,7) = 32.0/45.0;
    refNjNk(7,8) = 4.0/9.0;
    refNjNk(8,1) = (-2.0)/15.0;
    refNjNk(8,2) = (-8.0)/45.0;
    refNjNk(8,3) = (-8.0)/45.0;
    refNjNk(8,4) = (-2.0)/15.0;
    refNjNk(8,5) = 4.0/9.0;
    refNjNk(8,6) = 16.0/45.0;
    refNjNk(8,7) = 4.0/9.0;
    refNjNk(8,8) = 32.0/45.0;

// scale to bring this into physical space
    refNjNk *= 0.5*dx*0.5*dy;

// face mass matrix (x-lower)
    refFaceNjNk_xl = Lucee::Matrix<double>(faceShape, start);
    refFaceNjNk_xl(1,1) = 4.0/15.0;
    refFaceNjNk_xl(1,2) = 2.0/15.0;
    refFaceNjNk_xl(1,3) = (-1.0)/15.0;
    refFaceNjNk_xl(2,1) = 0;
    refFaceNjNk_xl(2,2) = 0;
    refFaceNjNk_xl(2,3) = 0;
    refFaceNjNk_xl(3,1) = 0;
    refFaceNjNk_xl(3,2) = 0;;
    refFaceNjNk_xl(3,3) = 0;
    refFaceNjNk_xl(4,1) = (-1.0)/15.0;
    refFaceNjNk_xl(4,2) = 2.0/15.0;
    refFaceNjNk_xl(4,3) = 4.0/15.0;;
    refFaceNjNk_xl(5,1) = 0;
    refFaceNjNk_xl(5,2) = 0;
    refFaceNjNk_xl(5,3) = 0;
    refFaceNjNk_xl(6,1) = 0;
    refFaceNjNk_xl(6,2) = 0;
    refFaceNjNk_xl(6,3) = 0;
    refFaceNjNk_xl(7,1) = 0;
    refFaceNjNk_xl(7,2) = 0;
    refFaceNjNk_xl(7,3) = 0;
    refFaceNjNk_xl(8,1) = 2.0/15.0;
    refFaceNjNk_xl(8,2) = 16.0/15.0;
    refFaceNjNk_xl(8,3) = 2.0/15.0;

// scale to bring this into physical space
    refFaceNjNk_xl *= 0.5*dy;

// face mass matrix (x-upper)
    refFaceNjNk_xu = Lucee::Matrix<double>(faceShape, start);
    refFaceNjNk_xu(1,1) = 0;
    refFaceNjNk_xu(1,2) = 0;
    refFaceNjNk_xu(1,3) = 0;
    refFaceNjNk_xu(2,1) = 4.0/15.0;
    refFaceNjNk_xu(2,2) = 2.0/15.0;
    refFaceNjNk_xu(2,3) = (-1.0)/15.0;
    refFaceNjNk_xu(3,1) = (-1.0)/15.0;
    refFaceNjNk_xu(3,2) = 2.0/15.0;
    refFaceNjNk_xu(3,3) = 4.0/15.0;
    refFaceNjNk_xu(4,1) = 0;
    refFaceNjNk_xu(4,2) = 0;
    refFaceNjNk_xu(4,3) = 0;
    refFaceNjNk_xu(5,1) = 0;
    refFaceNjNk_xu(5,2) = 0;;
    refFaceNjNk_xu(5,3) = 0;
    refFaceNjNk_xu(6,1) = 2.0/15.0;
    refFaceNjNk_xu(6,2) = 16.0/15.0;
    refFaceNjNk_xu(6,3) = 2.0/15.0;
    refFaceNjNk_xu(7,1) = 0;
    refFaceNjNk_xu(7,2) = 0;
    refFaceNjNk_xu(7,3) = 0;
    refFaceNjNk_xu(8,1) = 0;
    refFaceNjNk_xu(8,2) = 0;
    refFaceNjNk_xu(8,3) = 0;

// scale to bring this into physical space
    refFaceNjNk_xu *= 0.5*dy;

// face mass matrix (y-lower)
    refFaceNjNk_yl = Lucee::Matrix<double>(faceShape, start);
    refFaceNjNk_yl(1,1) = 4.0/15.0;
    refFaceNjNk_yl(1,2) = 2.0/15.0;
    refFaceNjNk_yl(1,3) = (-1.0)/15.0;
    refFaceNjNk_yl(2,1) = (-1.0)/15.0;
    refFaceNjNk_yl(2,2) = 2.0/15.0;
    refFaceNjNk_yl(2,3) = 4.0/15.0;
    refFaceNjNk_yl(3,1) = 0;
    refFaceNjNk_yl(3,2) = 0;
    refFaceNjNk_yl(3,3) = 0;
    refFaceNjNk_yl(4,1) = 0;
    refFaceNjNk_yl(4,2) = 0;
    refFaceNjNk_yl(4,3) = 0;
    refFaceNjNk_yl(5,1) = 2.0/15.0;
    refFaceNjNk_yl(5,2) = 16.0/15.0;
    refFaceNjNk_yl(5,3) = 2.0/15.0;
    refFaceNjNk_yl(6,1) = 0;
    refFaceNjNk_yl(6,2) = 0;
    refFaceNjNk_yl(6,3) = 0;
    refFaceNjNk_yl(7,1) = 0;
    refFaceNjNk_yl(7,2) = 0;
    refFaceNjNk_yl(7,3) = 0;
    refFaceNjNk_yl(8,1) = 0;
    refFaceNjNk_yl(8,2) = 0;
    refFaceNjNk_yl(8,3) = 0;

// scale to bring this into physical space
    refFaceNjNk_yl *= 0.5*dx;

// face mass matrix (y-upper)
    refFaceNjNk_yu = Lucee::Matrix<double>(faceShape, start);
    refFaceNjNk_yu(1,1) = 0;
    refFaceNjNk_yu(1,2) = 0;
    refFaceNjNk_yu(1,3) = 0;
    refFaceNjNk_yu(2,1) = 0;;
    refFaceNjNk_yu(2,2) = 0;
    refFaceNjNk_yu(2,3) = 0;
    refFaceNjNk_yu(3,1) = (-1.0)/15.0;
    refFaceNjNk_yu(3,2) = 2.0/15.0;
    refFaceNjNk_yu(3,3) = 4.0/15.0;
    refFaceNjNk_yu(4,1) = 4.0/15.0;
    refFaceNjNk_yu(4,2) = 2.0/15.0;;
    refFaceNjNk_yu(4,3) = (-1.0)/15.0;
    refFaceNjNk_yu(5,1) = 0;
    refFaceNjNk_yu(5,2) = 0;
    refFaceNjNk_yu(5,3) = 0;
    refFaceNjNk_yu(6,1) = 0;
    refFaceNjNk_yu(6,2) = 0;
    refFaceNjNk_yu(6,3) = 0;
    refFaceNjNk_yu(7,1) = 2.0/15.0;
    refFaceNjNk_yu(7,2) = 16.0/15.0;
    refFaceNjNk_yu(7,3) = 2.0/15.0;
    refFaceNjNk_yu(8,1) = 0;
    refFaceNjNk_yu(8,2) = 0;
    refFaceNjNk_yu(8,3) = 0;

// scale to bring this into physical space
    refFaceNjNk_yu *= 0.5*dx;

// stiffness matrix (automatically generated. See scripts/serendipity-2D.mac)
    refDNjDNk = Lucee::Matrix<double>(shape, start);
    refDNjDNk(1,1) = (373*dy2+328*dx2)/(dx2*dy2)/180.0+(43*dy2+88*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(1,2) = -(-(187*dy2+68*dx2)/6.0-(37*dy2+68*dx2)/6.0)/(dx2*dy2)/30.0;
    refDNjDNk(1,3) = 2.0*(23*dy2+23*dx2)/(45.0*dx2*dy2);
    refDNjDNk(1,4) = 2.0*(17*dy2+28*dx2)/(45.0*dx2*dy2);
    refDNjDNk(1,5) = -(140*dy2-36*dx2)/(dx2*dy2)/45.0-(20*dy2+24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(1,6) = (21*dy2-160*dx2)/(dx2*dy2)/180.0-(69*dy2+160*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(1,7) = -(40*dy2+36*dx2)/(dx2*dy2)/45.0-(40*dy2-24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(1,8) = (69*dy2-440*dx2)/(dx2*dy2)/180.0-(21*dy2+200*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(2,1) = -(-(187*dy2+68*dx2)/6.0-(37*dy2+68*dx2)/6.0)/(dx2*dy2)/30.0;
    refDNjDNk(2,2) = (373*dy2+328*dx2)/(dx2*dy2)/180.0+(43*dy2+88*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(2,3) = 2.0*(17*dy2+28*dx2)/(45.0*dx2*dy2);
    refDNjDNk(2,4) = 2.0*(23*dy2+23*dx2)/(45.0*dx2*dy2);
    refDNjDNk(2,5) = -(140*dy2-36*dx2)/(dx2*dy2)/45.0-(20*dy2+24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(2,6) = (69*dy2-440*dx2)/(dx2*dy2)/180.0-(21*dy2+200*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(2,7) = -(40*dy2+36*dx2)/(dx2*dy2)/45.0-(40*dy2-24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(2,8) = (21*dy2-160*dx2)/(dx2*dy2)/180.0-(69*dy2+160*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(3,1) = 2.0*(23*dy2+23*dx2)/(45.0*dx2*dy2);
    refDNjDNk(3,2) = 2.0*(17*dy2+28*dx2)/(45.0*dx2*dy2);
    refDNjDNk(3,3) = (373*dy2+328*dx2)/(dx2*dy2)/180.0+(43*dy2+88*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(3,4) = -(-(187*dy2+68*dx2)/6.0-(37*dy2+68*dx2)/6.0)/(dx2*dy2)/30.0;
    refDNjDNk(3,5) = -(40*dy2+36*dx2)/(dx2*dy2)/45.0-(40*dy2-24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(3,6) = (69*dy2-440*dx2)/(dx2*dy2)/180.0-(21*dy2+200*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(3,7) = -(140*dy2-36*dx2)/(dx2*dy2)/45.0-(20*dy2+24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(3,8) = (21*dy2-160*dx2)/(dx2*dy2)/180.0-(69*dy2+160*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(4,1) = 2.0*(17*dy2+28*dx2)/(45.0*dx2*dy2);
    refDNjDNk(4,2) = 2.0*(23*dy2+23*dx2)/(45.0*dx2*dy2);
    refDNjDNk(4,3) = -(-(187*dy2+68*dx2)/6.0-(37*dy2+68*dx2)/6.0)/(dx2*dy2)/30.0;
    refDNjDNk(4,4) = (373*dy2+328*dx2)/(dx2*dy2)/180.0+(43*dy2+88*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(4,5) = -(40*dy2+36*dx2)/(dx2*dy2)/45.0-(40*dy2-24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(4,6) = (21*dy2-160*dx2)/(dx2*dy2)/180.0-(69*dy2+160*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(4,7) = -(140*dy2-36*dx2)/(dx2*dy2)/45.0-(20*dy2+24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(4,8) = (69*dy2-440*dx2)/(dx2*dy2)/180.0-(21*dy2+200*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(5,1) = -(140*dy2-36*dx2)/(dx2*dy2)/45.0-(20*dy2+24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(5,2) = -(140*dy2-36*dx2)/(dx2*dy2)/45.0-(20*dy2+24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(5,3) = -(40*dy2+36*dx2)/(dx2*dy2)/45.0-(40*dy2-24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(5,4) = -(40*dy2+36*dx2)/(dx2*dy2)/45.0-(40*dy2-24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(5,5) = 2.0*((140*dy2+24*dx2)/3.0+(20*dy2+24*dx2)/3.0)/(15.0*dx2*dy2);
    refDNjDNk(5,6) = 0;
    refDNjDNk(5,7) = 4.0*(40*dy2-24*dx2)/(45.0*dx2*dy2);
    refDNjDNk(5,8) = 0;
    refDNjDNk(6,1) = (21*dy2-160*dx2)/(dx2*dy2)/180.0-(69*dy2+160*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(6,2) = (69*dy2-440*dx2)/(dx2*dy2)/180.0-(21*dy2+200*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(6,3) = (69*dy2-440*dx2)/(dx2*dy2)/180.0-(21*dy2+200*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(6,4) = (21*dy2-160*dx2)/(dx2*dy2)/180.0-(69*dy2+160*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(6,5) = 0;
    refDNjDNk(6,6) = 2.0*(48*dy2+160*dx2)/(45.0*dx2*dy2);
    refDNjDNk(6,7) = 0;
    refDNjDNk(6,8) = (-4.0)*(24*dy2-40*dx2)/(45.0*dx2*dy2);
    refDNjDNk(7,1) = -(40*dy2+36*dx2)/(dx2*dy2)/45.0-(40*dy2-24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(7,2) = -(40*dy2+36*dx2)/(dx2*dy2)/45.0-(40*dy2-24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(7,3) = -(140*dy2-36*dx2)/(dx2*dy2)/45.0-(20*dy2+24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(7,4) = -(140*dy2-36*dx2)/(dx2*dy2)/45.0-(20*dy2+24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(7,5) = 4.0*(40*dy2-24*dx2)/(45.0*dx2*dy2);
    refDNjDNk(7,6) = 0;
    refDNjDNk(7,7) = 2.0*((140*dy2+24*dx2)/3.0+(20*dy2+24*dx2)/3.0)/(15.0*dx2*dy2);
    refDNjDNk(7,8) = 0;
    refDNjDNk(8,1) = (69*dy2-440*dx2)/(dx2*dy2)/180.0-(21*dy2+200*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(8,2) = (21*dy2-160*dx2)/(dx2*dy2)/180.0-(69*dy2+160*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(8,3) = (21*dy2-160*dx2)/(dx2*dy2)/180.0-(69*dy2+160*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(8,4) = (69*dy2-440*dx2)/(dx2*dy2)/180.0-(21*dy2+200*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(8,5) = 0;
    refDNjDNk(8,6) = (-4.0)*(24*dy2-40*dx2)/(45.0*dx2*dy2);
    refDNjDNk(8,7) = 0;
    refDNjDNk(8,8) = 2.0*(48*dy2+160*dx2)/(45.0*dx2*dy2);

// scale to bring this into physical space
    refDNjDNk *= 0.5*dx*0.5*dy;

// grad-stiffness matrix (automatically generated. See scripts/serendipity-2D.mac)
    refDNjNk_0 = Lucee::Matrix<double>(shape, start);
    refDNjNk_0(1,1) = (-4.0)/(15.0*dx);
    refDNjNk_0(1,2) = 8.0/(45.0*dx);
    refDNjNk_0(1,3) = 1/dx/15.0;
    refDNjNk_0(1,4) = 1/dx/15.0;
    refDNjNk_0(1,5) = (-4.0)/(9.0*dx);
    refDNjNk_0(1,6) = 14.0/(45.0*dx);
    refDNjNk_0(1,7) = 0;
    refDNjNk_0(1,8) = (-26.0)/(45.0*dx);
    refDNjNk_0(2,1) = (-8.0)/(45.0*dx);
    refDNjNk_0(2,2) = 4.0/(15.0*dx);
    refDNjNk_0(2,3) = -1/dx/15.0;
    refDNjNk_0(2,4) = -1/dx/15.0;
    refDNjNk_0(2,5) = 4.0/(9.0*dx);
    refDNjNk_0(2,6) = 26.0/(45.0*dx);
    refDNjNk_0(2,7) = 0;
    refDNjNk_0(2,8) = (-14.0)/(45.0*dx);
    refDNjNk_0(3,1) = -1/dx/15.0;
    refDNjNk_0(3,2) = -1/dx/15.0;
    refDNjNk_0(3,3) = 4.0/(15.0*dx);
    refDNjNk_0(3,4) = (-8.0)/(45.0*dx);
    refDNjNk_0(3,5) = 0;
    refDNjNk_0(3,6) = 26.0/(45.0*dx);
    refDNjNk_0(3,7) = 4.0/(9.0*dx);
    refDNjNk_0(3,8) = (-14.0)/(45.0*dx);
    refDNjNk_0(4,1) = 1/dx/15.0;
    refDNjNk_0(4,2) = 1/dx/15.0;
    refDNjNk_0(4,3) = 8.0/(45.0*dx);
    refDNjNk_0(4,4) = (-4.0)/(15.0*dx);
    refDNjNk_0(4,5) = 0;
    refDNjNk_0(4,6) = 14.0/(45.0*dx);
    refDNjNk_0(4,7) = (-4.0)/(9.0*dx);
    refDNjNk_0(4,8) = (-26.0)/(45.0*dx);
    refDNjNk_0(5,1) = 4.0/(9.0*dx);
    refDNjNk_0(5,2) = (-4.0)/(9.0*dx);
    refDNjNk_0(5,3) = 0;
    refDNjNk_0(5,4) = 0;
    refDNjNk_0(5,5) = 0;
    refDNjNk_0(5,6) = (-8.0)/(9.0*dx);
    refDNjNk_0(5,7) = 0;
    refDNjNk_0(5,8) = 8.0/(9.0*dx);
    refDNjNk_0(6,1) = (-14.0)/(45.0*dx);
    refDNjNk_0(6,2) = (-14.0)/(45.0*dx);
    refDNjNk_0(6,3) = (-14.0)/(45.0*dx);
    refDNjNk_0(6,4) = (-14.0)/(45.0*dx);
    refDNjNk_0(6,5) = 8.0/(9.0*dx);
    refDNjNk_0(6,6) = 16.0/(15.0*dx);
    refDNjNk_0(6,7) = 8.0/(9.0*dx);
    refDNjNk_0(6,8) = 16.0/(15.0*dx);
    refDNjNk_0(7,1) = 0;
    refDNjNk_0(7,2) = 0;
    refDNjNk_0(7,3) = (-4.0)/(9.0*dx);
    refDNjNk_0(7,4) = 4.0/(9.0*dx);
    refDNjNk_0(7,5) = 0;
    refDNjNk_0(7,6) = (-8.0)/(9.0*dx);
    refDNjNk_0(7,7) = 0;
    refDNjNk_0(7,8) = 8.0/(9.0*dx);
    refDNjNk_0(8,1) = 14.0/(45.0*dx);
    refDNjNk_0(8,2) = 14.0/(45.0*dx);
    refDNjNk_0(8,3) = 14.0/(45.0*dx);
    refDNjNk_0(8,4) = 14.0/(45.0*dx);
    refDNjNk_0(8,5) = (-8.0)/(9.0*dx);
    refDNjNk_0(8,6) = (-16.0)/(15.0*dx);
    refDNjNk_0(8,7) = (-8.0)/(9.0*dx);
    refDNjNk_0(8,8) = (-16.0)/(15.0*dx);

// scale to bring this into physical space
    refDNjNk_0 *= 0.5*dx*0.5*dy;

// grad-stiffness matrix (automatically generated. See scripts/serendipity-2D.mac)
    refDNjNk_1 = Lucee::Matrix<double>(shape, start);
    refDNjNk_1(1,1) = (-4.0)/(15.0*dy);
    refDNjNk_1(1,2) = 1/dy/15.0;
    refDNjNk_1(1,3) = 1/dy/15.0;
    refDNjNk_1(1,4) = 8.0/(45.0*dy);
    refDNjNk_1(1,5) = (-26.0)/(45.0*dy);
    refDNjNk_1(1,6) = 0;
    refDNjNk_1(1,7) = 14.0/(45.0*dy);
    refDNjNk_1(1,8) = (-4.0)/(9.0*dy);
    refDNjNk_1(2,1) = 1/dy/15.0;
    refDNjNk_1(2,2) = (-4.0)/(15.0*dy);
    refDNjNk_1(2,3) = 8.0/(45.0*dy);
    refDNjNk_1(2,4) = 1/dy/15.0;
    refDNjNk_1(2,5) = (-26.0)/(45.0*dy);
    refDNjNk_1(2,6) = (-4.0)/(9.0*dy);
    refDNjNk_1(2,7) = 14.0/(45.0*dy);
    refDNjNk_1(2,8) = 0;
    refDNjNk_1(3,1) = -1/dy/15.0;
    refDNjNk_1(3,2) = (-8.0)/(45.0*dy);
    refDNjNk_1(3,3) = 4.0/(15.0*dy);
    refDNjNk_1(3,4) = -1/dy/15.0;
    refDNjNk_1(3,5) = (-14.0)/(45.0*dy);
    refDNjNk_1(3,6) = 4.0/(9.0*dy);
    refDNjNk_1(3,7) = 26.0/(45.0*dy);
    refDNjNk_1(3,8) = 0;
    refDNjNk_1(4,1) = (-8.0)/(45.0*dy);
    refDNjNk_1(4,2) = -1/dy/15.0;
    refDNjNk_1(4,3) = -1/dy/15.0;
    refDNjNk_1(4,4) = 4.0/(15.0*dy);
    refDNjNk_1(4,5) = (-14.0)/(45.0*dy);
    refDNjNk_1(4,6) = 0;
    refDNjNk_1(4,7) = 26.0/(45.0*dy);
    refDNjNk_1(4,8) = 4.0/(9.0*dy);
    refDNjNk_1(5,1) = 14.0/(45.0*dy);
    refDNjNk_1(5,2) = 14.0/(45.0*dy);
    refDNjNk_1(5,3) = 14.0/(45.0*dy);
    refDNjNk_1(5,4) = 14.0/(45.0*dy);
    refDNjNk_1(5,5) = (-16.0)/(15.0*dy);
    refDNjNk_1(5,6) = (-8.0)/(9.0*dy);
    refDNjNk_1(5,7) = (-16.0)/(15.0*dy);
    refDNjNk_1(5,8) = (-8.0)/(9.0*dy);
    refDNjNk_1(6,1) = 0;
    refDNjNk_1(6,2) = 4.0/(9.0*dy);
    refDNjNk_1(6,3) = (-4.0)/(9.0*dy);
    refDNjNk_1(6,4) = 0;
    refDNjNk_1(6,5) = 8.0/(9.0*dy);
    refDNjNk_1(6,6) = 0;
    refDNjNk_1(6,7) = (-8.0)/(9.0*dy);
    refDNjNk_1(6,8) = 0;
    refDNjNk_1(7,1) = (-14.0)/(45.0*dy);
    refDNjNk_1(7,2) = (-14.0)/(45.0*dy);
    refDNjNk_1(7,3) = (-14.0)/(45.0*dy);
    refDNjNk_1(7,4) = (-14.0)/(45.0*dy);
    refDNjNk_1(7,5) = 16.0/(15.0*dy);
    refDNjNk_1(7,6) = 8.0/(9.0*dy);
    refDNjNk_1(7,7) = 16.0/(15.0*dy);
    refDNjNk_1(7,8) = 8.0/(9.0*dy);
    refDNjNk_1(8,1) = 4.0/(9.0*dy);
    refDNjNk_1(8,2) = 0;
    refDNjNk_1(8,3) = 0;
    refDNjNk_1(8,4) = (-4.0)/(9.0*dy);
    refDNjNk_1(8,5) = 8.0/(9.0*dy);
    refDNjNk_1(8,6) = 0;
    refDNjNk_1(8,7) = (-8.0)/(9.0*dy);
    refDNjNk_1(8,8) = 0;

// scale to bring this into physical space
    refDNjNk_1 *= 0.5*dx*0.5*dy;
  }

  void
  SerendipityElement2D::getGlobalIndices(int ix, int iy, std::vector<int>& glob,
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

    if ((ix<numX) && (iy<numY))
    { // nodes inside grid proper
      glob.resize(3);
      glob[0] = F_func(numX, numY, ix, iy); // node 1
      glob[1] = F_func(numX, numY, ix, iy) + 1; // node 5
      glob[2] = G_func(numX, numY, ix, iy); // node 8

      loc.resize(3);
      loc[0] = 0; loc[1] = 1; loc[2] = 2;
    }
    if ((ix==numX) && (iy<numY))
    { // right edge nodes
      glob.resize(2);
      glob[0] = F_func(numX, numY, ix, iy); // node 1
      glob[1] = G_func(numX, numY, ix, iy); // node 8

      loc.resize(2);
      loc[0] = 0; loc[1] = 2;
    }
    if ((ix<numX) && (iy==numY))
    { // top edge nodes
      glob.resize(2);
      glob[0] = F_func(numX, numY, ix, iy); // node 1
      glob[1] = F_func(numX, numY, ix, iy) + 1; // node 5

      loc.resize(2);
      loc[0] = 0; loc[1] = 1;
    }
    if ((ix==numX) && (iy==numY))
    { // top right corner
      glob.resize(1);
      glob[0] = F_func(numX, numY, ix, iy); // node 1

      loc.resize(1);
      loc[0] = 0;
    }
  }
}
