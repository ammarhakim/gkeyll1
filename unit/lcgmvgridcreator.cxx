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
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>

// lucee includes
#include <LcCmdLineArgs.h>
#include <LcGmvGridCreator.h>
#include <LcTest.h>
#include <LcUnstructGrid.h>

void
test_1(const Lucee::UnstructGrid<double>& grid)
{
  double xv[3];
// create iterator over nodes
  Lucee::UnstructGrid<double>::ElemIterator<0> vitr(grid);
  for ( ; !vitr.atEnd(); ++vitr)
  {
// get nodal coordinates
    vitr->fillWithCoordinates(xv);
  }

// create iterator over cells
  Lucee::UnstructGrid<double>::ElemIterator<3> citr(grid);

// create incidence iterator
  Lucee::UnstructGrid<double>::IncidenceIterator<3, 0> c2vItr(grid);
}

int
main (int argc, char *argv[])
{
  LC_BEGIN_TESTS("lcgmvgridcreator");

  Lucee::CmdLineArgs cmd("lcgmvgridcreator");
  cmd.addArg("i", "GMVFILE", "GMV file to read");
  cmd.addArg("d", "NDIM", "Grid dimension");
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
  unsigned ndim = 3;
  if (cmd.hasArg("d"))
    ndim = std::atoi(cmd.getArg("d").c_str());

// open file
  std::ifstream gmvFile(gmvFileNm.c_str());
  if (!gmvFile)
  {
    std::cerr << "Unable to open GMV file " << gmvFileNm << std::endl;
    exit(1);
  }
// create from reader
  Lucee::GmvGridCreator<double> gmvRdr(ndim, gmvFile);
// create a new unstructured grid
  Lucee::UnstructGrid<double> ugrid;
  ugrid.constructFromCreator(gmvRdr);

// write grid to HDF5 file
  ugrid.write("ugrid.h5");

// ensure correct number of nodes and cells
  LC_ASSERT("Checking number of nodes", ugrid.getNumVertices() == 7856);
  LC_ASSERT("Checking number of cells", ugrid.getNumCells() == 39597);

// now run other tests
  test_1(ugrid);

  LC_END_TESTS;
}
