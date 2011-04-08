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

// lucee includes
#include <LcUnstructGrid.h>

namespace Lucee
{
// set ids for grid creators
  template <> const char *UnstructGrid<double>::id = "Unstruct";

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
  UnstructGrid<REAL>::getNumTri() const
  {
    std::map<short, unsigned>::const_iterator itr
      = cellCount.find(TRI_CELL_T);
    return itr->second;
  }

  template <typename REAL>
  unsigned
  UnstructGrid<REAL>::getNumQuad() const
  {
    std::map<short, unsigned>::const_iterator itr
      = cellCount.find(QUAD_CELL_T);
    return itr->second;
  }

  template <typename REAL>
  unsigned
  UnstructGrid<REAL>::getNumTet() const
  {
    std::map<short, unsigned>::const_iterator itr
      = cellCount.find(TET_CELL_T);
    return itr->second;
  }

  template <typename REAL>
  unsigned
  UnstructGrid<REAL>::getNumHex() const
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

    unsigned nn; // for 2D
// write out connectivity
    if (ndim == 2)
      nn = 3+1; // tri
    else
      nn = 4+1; // tet
    std::vector<int> conn(getNumCells()*nn);
    
// now fill up connectivities
    std::vector<Lucee::UnstructConnectivity>::const_iterator c2vItr
      = connectivity.begin()+4*ndim; // iterator to cell->0 connectivity
    for (unsigned ic=0; ic<getNumCells(); ++ic)
    {
      conn[nn*ic+0] = nn-1;
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

    ctor.fillWithConnectivity(connectivity[4*ndim+0]);
    ddprime[4*ndim+0] = true; // set flag as now ndim->0 connectivity is stored

    cellCount[TRI_CELL_T] = ctor.getNumTri();
    cellCount[QUAD_CELL_T] = ctor.getNumQuad();
    cellCount[TET_CELL_T] = ctor.getNumTet();
    cellCount[HEX_CELL_T] = ctor.getNumHex();
  }

// instantiations
  template class UnstructGrid<float>;
  template class UnstructGrid<double>;
}
