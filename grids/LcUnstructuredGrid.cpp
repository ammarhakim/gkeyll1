/*
 * LcUnstructuredGrid.cpp
 *
 *  Created on: May 18, 2015
 *      Author: john
 */

#include "LcUnstructuredGrid.h"

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCartProdDecompRegionCalc.h>
#include <LcGlobals.h>
#include <LcPointerHolder.h>
#include <LcUnstructuredGrid.h>

// txbase includes
#include <TxCommBase.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  template<> const char *UnstructuredGrid<2>::id = "Unstructured2D";
  template<> const char *UnstructuredGrid<3>::id = "Unstructured3D";

  template<unsigned NDIM>
  UnstructuredGrid<NDIM>::~UnstructuredGrid()
  {
  }

  template<unsigned NDIM>
  void UnstructuredGrid<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// read number of cells in global region
    std::string filename = tbl.getString("filename");

    mb = new (std::nothrow) moab::Core;

    moab::ErrorCode rval = mb->load_mesh(filename.c_str());

    mb->get_entities_by_dimension(0, NDIM - 1, faces);
    mb->get_entities_by_dimension(0, NDIM, cells);
    mb->get_entities_by_type(0, moab::MBVERTEX, vertices);

// get comm pointers
    TxCommBase * comm = Loki::SingletonHolder<Lucee::Globals>::Instance().comm;
    TxCommBase * momComm =
        Loki::SingletonHolder<Lucee::Globals>::Instance().comm;

// set valid communicators for grid
    this->setComm(comm);
    this->setMomComm(momComm);
// set I/O flag for safe ranks
    this->setIsSafeToWrite(true);
// compute local region
    //localRgn = decompRgn->getRegion(comm->getRank());
  }

  template<unsigned NDIM>
  TxIoNodeType UnstructuredGrid<NDIM>::writeToFile(TxIoBase& io, TxIoNodeType& node,
      const std::string& nm) {
    return node;
  }

  template<unsigned NDIM>
  unsigned UnstructuredGrid<NDIM>::getNumCells(unsigned dir) const
  {
    return cells.size();
  }

  template<unsigned NDIM>
  void UnstructuredGrid<NDIM>::setIndex(int i) const
  {
    currIdx[0] = i;
  }

  template<unsigned NDIM>
  void UnstructuredGrid<NDIM>::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
    //lfm.appendFunc("filename", luaGetFilename);
  }

  template<unsigned NDIM>
  UnstructuredGrid<NDIM>::UnstructuredGrid()
  {
    mb = NULL;
  }

  template<unsigned NDIM>
  UnstructuredGrid<NDIM>&
  UnstructuredGrid<NDIM>::operator=(const UnstructuredGrid<NDIM>& sg)
  {
    throw Lucee::Except("UnstructuredGrid::operator=: Not implemented");
    return *this;
  }

// instantiations

  template class UnstructuredGrid<2> ;
  template class UnstructuredGrid<3> ;
}
