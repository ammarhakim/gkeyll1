/**
 * @file	lchdf5io.cxx
 *
 * @brief	Unit tests for Lucee::Hdf5Io class
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// std includes
#include <iostream>
#include <cstdlib>
#include <fstream>

// lucee includes
#include <LcCmdLineArgs.h>
#include <LcGmvGridCreator.h>
#include <LcTest.h>
#include <LcUnstructGrid.h>

int
main (int argc, char *argv[])
{
  LC_BEGIN_TESTS("lcgmvgridcreator");

  Lucee::CmdLineArgs cmd("lcgmvgridcreator");
  cmd.addArg("i", "GMVFILE", "GMV file to read");
// parse it
  cmd.parse(argc, argv);
// show help if requested
  if (cmd.hasSwitch("h"))
  {
    cmd.showHelp();
    exit(1);
  }

  std::string gmvFileNm;
  if (cmd.hasArg("i"))
    gmvFileNm = cmd.getArg("i");
  else
  {
    std::cerr << "Must provide a GMV input file to read!" << std::endl;
    cmd.showHelp();
    exit(1);
  }

// open file
  std::ifstream gmvFile(gmvFileNm.c_str());
  if (!gmvFile)
  {
    std::cerr << "Unable to open GMV file " << gmvFileNm << std::endl;
    exit(1);
  }
// create from reader
  Lucee::GmvGridCreator<double> gmvRdr(3, gmvFile);
// create a new unstructured grid
  Lucee::UnstructGrid<double> ugrid;
  ugrid.constructFromCreator(gmvRdr);

  double xv[3];
// create an iterator over vertices
  Lucee::UnstructGrid<double>::ElemIterator<0> itr(ugrid);
  for ( ; !itr.atEnd(); ++itr)
  {
    itr->fillWithCoordinates(xv);
    for (unsigned i=0; i<3; ++i)
      std::cout << xv[i] << " ";
    std::cout << std::endl;
  }

// write grid to HDF5 file
  ugrid.write("ugrid.h5");

// ensure correct number of nodes and cells
  LC_ASSERT("Checking number of nodes", ugrid.getNumVertices() == 7856);
  LC_ASSERT("Checking number of cells", ugrid.getNumCells() == 39597);

  LC_END_TESTS;
}
