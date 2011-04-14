/**
 * @file	LcUnstructGrid.cpp
 *
 * @brief	Unstructured grid class.
 *
 * @version	$Id$
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// std includes
#include <iostream>
#include <valarray>
#include <vector>

// lucee includes
#include <LcUnstructGrid.h>

namespace Lucee
{
// set ids for grid creators
  template <> const char *UnstructGrid<double>::id = "Unstruct";

/**
 * Convenient method to map connectivity indices to linear index.
 *
 * @param d Dimension to connect (d->dprime).
 * @param dprime Dimension to connect to (d->dprime).
 * @return linear index.
 */
  inline static unsigned
  getConnIndex(unsigned d, unsigned dprime)
  {
    return 4*d+dprime;
  }

  template <typename REAL>
  UnstructGrid<REAL>::UnstructGrid()
    : ndim(3), ddprime(4*4), connectivity(4*4)
  {
    this->setName("Unstruct");
// initialize ddprime as no connections have been made
    for (unsigned i=0; i<16; ++i)
      ddprime[i] = false;
  }

  template <typename REAL>
  void
  UnstructGrid<REAL>::readInput(Lucee::LuaTable& tbl)
  {
    
  }

  template <typename REAL>
  unsigned
  UnstructGrid<REAL>::getNumVertices() const
  {
    return geometry.vcoords.size()/3;
  }

  template <typename REAL>
  unsigned
  UnstructGrid<REAL>::getNumCells() const
  {
    return connectivity[4*ndim+0].offsets.size()-1;
  }

  template <typename REAL>
  unsigned
  UnstructGrid<REAL>::getNumTriCells() const
  {
    std::map<short, unsigned>::const_iterator itr
      = cellCount.find(TRI_CELL_T);
    return itr->second;
  }

  template <typename REAL>
  unsigned
  UnstructGrid<REAL>::getNumQuadCells() const
  {
    std::map<short, unsigned>::const_iterator itr
      = cellCount.find(QUAD_CELL_T);
    return itr->second;
  }

  template <typename REAL>
  unsigned
  UnstructGrid<REAL>::getNumTetCells() const
  {
    std::map<short, unsigned>::const_iterator itr
      = cellCount.find(TET_CELL_T);
    return itr->second;
  }

  template <typename REAL>
  unsigned
  UnstructGrid<REAL>::getNumHexCells() const
  {
    std::map<short, unsigned>::const_iterator itr
      = cellCount.find(HEX_CELL_T);
    return itr->second;
  }

  template <typename REAL>
  Lucee::IoNodeType
  UnstructGrid<REAL>::writeToFile(Lucee::IoBase& io, Lucee::IoNodeType& node,
    const std::string& nm)
  {
// create a group to write grid data
    Lucee::IoNodeType gn = io.createGroup(node, nm);
// add VizSchema markup
    io.writeStrAttribute(gn, "vsType", "mesh");
    io.writeStrAttribute(gn, "vsKind", "unstructured");
    io.writeStrAttribute(gn, "vsPoints", "vertices");
    io.writeStrAttribute(gn, "vsPolygons", "cells");

    std::vector<size_t> dsSize(2);
    std::vector<size_t> dsBeg(2);
    std::vector<size_t> dsLen(2);
    
// write vertices to file
    dsLen[0] = dsSize[0] = getNumVertices();
    dsLen[1] = dsSize[1] = 3;
    io.writeDataSet(gn, "vertices", dsSize, dsBeg, dsLen, &geometry.vcoords[0]);

    unsigned nn;
// determine maximum node count in cell
    if (ndim==1)
    {
      throw Lucee::Except("UnstructGrid::writeDataSet: 1D unstructured mesh not currently supported");
    }
    else if (ndim==2)
    { // MUST preserve order of following ifs (in increasing order of node count)
      if (getNumTriCells() > 0) nn = 3;
      if (getNumQuadCells() > 0) nn = 4;
    }
    else
    { // MUST preserve order of following ifs (in increasing order of node count)
      if (getNumTetCells() > 0) nn = 4;
      if (getNumHexCells() > 0) nn = 8;
    }
    nn = nn+1; // one extra to write number of connected nodes per cell
// allocate array to write cell connectivities
    std::valarray<int> conn(getNumCells()*nn);
    conn = 0;
    
// fill connectivities
    std::vector<Lucee::UnstructConnectivity>::const_iterator c2vItr
      = connectivity.begin()+getConnIndex(ndim,0); // set itr to cell->0 connectivity
    for (unsigned ic=0; ic<getNumCells(); ++ic)
    {
      conn[nn*ic+0] = c2vItr->offsets[ic+1]-c2vItr->offsets[ic];
      unsigned count=1;
      for (unsigned iv=c2vItr->offsets[ic]; iv<c2vItr->offsets[ic+1]; ++iv)
        conn[nn*ic+count++] = c2vItr->indices[iv];
    }

// write connectivities
    dsLen[0] = dsSize[0] = getNumCells();
    dsLen[1] = dsSize[1] = nn;
    io.writeDataSet(gn, "cells", dsSize, dsBeg, dsLen, &conn[0]);

    return gn;
  }

  template <typename REAL>
  void
  UnstructGrid<REAL>::constructFromCreator(const Lucee::UnstructGridCreator<REAL>& ctor)
  {
    ndim = ctor.getDim();

    ctor.fillWithGeometry(geometry);
    ddprime[getConnIndex(ndim, 0)] = true; // set flag as ndim->0 connectivity is always stored
    ctor.fillWithConnectivity(connectivity[getConnIndex(ndim, 0)]);

// compute cell volume and centroids (TODO)
    geometry.setNumCells(getNumCells());

    cellCount[TRI_CELL_T] = ctor.getNumTri();
    cellCount[QUAD_CELL_T] = ctor.getNumQuad();
    cellCount[TET_CELL_T] = ctor.getNumTet();
    cellCount[HEX_CELL_T] = ctor.getNumHex();
  }

  template <typename REAL>
  const Lucee::UnstructConnectivity&
  UnstructGrid<REAL>::getConnectivity(unsigned d, unsigned dprime) const
  {
// check if connectivity exists
    if (ddprime[getConnIndex(d, dprime)])
      return connectivity[getConnIndex(d, dprime)];
    throw Lucee::Except("UnstructGrid::getConnectivity: Connectivity not computed!");
  }

// instantiations
  template class UnstructGrid<float>;
  template class UnstructGrid<double>;
}
