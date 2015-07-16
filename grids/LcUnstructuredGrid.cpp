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
#include <moab/FileOptions.hpp>

namespace Lucee
{
  template<> const char *UnstructuredGrid<1>::id = "Unstructured1D";
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
    filename = tbl.getString("filename");
    writeExtension = tbl.getString("writeExtension");

    //outFilename = tbl.getString("outFilename");

    mb = new (std::nothrow) moab::Core;

    mb->create_meshset(moab::MESHSET_SET, set);
    moab::ErrorCode rval = mb->load_file(filename.c_str(), &set, readOpts);

    mb->get_entities_by_dimension(0, NDIM - 1, faces);
    mb->get_entities_by_dimension(0, NDIM, cells);
    mb->get_entities_by_type(0, moab::MBVERTEX, vertices);

    std::cout << "numCells " << cells.size() << "\n";
    std::cout << "numFaces " << faces.size() << "\n";
    std::cout << "numVertices " << vertices.size() << "\n";

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
  void UnstructuredGrid<NDIM>::write(const std::string& nm)
  {
    // output prefix
    std::string outPrefix =
        Loki::SingletonHolder<Lucee::Globals>::Instance().outPrefix;
    std::string outNm = outPrefix + "_" + nm;

    //std::string newName = nm;
    outNm.append(".").append(writeExtension);

    std::cout << "newName is " << outNm << "\n";

    mb->write_file(outNm.c_str());

    /*if (writeExtension == "cgns")
    {
      mb->write_file(outNm.c_str(), "CGNS");
    } else if (writeExtension == "vtk")
    {
      mb->write_file(outNm.c_str(), "VTK");
    } else if (writeExtension == "gmv")
    {
      mb->write_file(outNm.c_str(), "GMV");
    } else if (writeExtension =="gmsh") {
      mb->write_file(outNm.c_str(), "GMSH");
    } else
    {
      std::cout << "No file was written!\n";
    }*/

  }

  //TODO: this actually isn't used now
  template<unsigned NDIM>
   TxIoNodeType UnstructuredGrid<NDIM>::writeToFile(TxIoBase& io,
   TxIoNodeType& node, const std::string& nm)
   {
   std::cout << "writing output\n";

   std::vector<std::basic_string<char>> dummy;


   //moab::ReaderWriterSet tSet(dynamic_cast<moab::Core*>(mb));
   //moab::WriterIface* wCGNS = tSet.get_file_extension_writer(nm);

   std::cout << "setting file options\n";

   //moab::FileOptions writeOpts(NULL);

   std::cout << "finally writing\n";

   std::cout << "writeExtension " << writeExtension << "\n";

   std::string newName = nm;
   newName.append(".").append(writeExtension);

   std::cout << "newName is " << newName << "\n";

   if(writeExtension=="cgns") {
   mb->write_file(newName.c_str(),"CGNS");
   } else if(writeExtension=="vtk") {
   mb->write_file(newName.c_str(),"VTK");
   } else if(writeExtension=="gmv") {
   mb->write_file(newName.c_str(),"GMV");
   } else {
   std::cout << "No file was written!\n";
   }

   //mb->write_file(nm.c_str());

   std::cout << "exiting this thing\n";

   //mb->write_file(nm.c_str(), 0, writeOpts, &set, 1);

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
    readOpts = NULL;
    //writeOpts = NULL;
  }

  template<unsigned NDIM>
  UnstructuredGrid<NDIM>&
  UnstructuredGrid<NDIM>::operator=(const UnstructuredGrid<NDIM>& sg)
  {
    throw Lucee::Except("UnstructuredGrid::operator=: Not implemented");
    return *this;
  }

// instantiations
  template class UnstructuredGrid<1> ;
  template class UnstructuredGrid<2> ;
  template class UnstructuredGrid<3> ;
}
