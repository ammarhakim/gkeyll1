/**
 * @file	LcSolverAssembly.cpp
 *
 * @brief	Solver that assembles updaters to create a simulation. *
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGenericFactory.h>
#include <LcHdf5Io.h>
#include <LcLogStream.h>
#include <LcLogger.h>
#include <LcObjCreator.h>
#include <LcSolverAssembly.h>

// std includes
#include <memory>

namespace Lucee
{
  const char *SolverAssembly::id = "Assembly";
  
  SolverAssembly::SolverAssembly()
    : SolverIfc(SolverAssembly::id)
  {
  }

  SolverAssembly::~SolverAssembly()
  {
// delete all grids
    std::map<std::string, Lucee::GridIfc*>::iterator gItr;
    for (gItr = gridMap.begin(); gItr != gridMap.end(); ++gItr)
      delete gItr->second;
    gridMap.clear();
// delete all data-structs
    std::map<std::string, Lucee::DataStructIfc*>::iterator dsItr;
    for (dsItr = dataStructMap.begin(); dsItr != dataStructMap.end(); ++dsItr)
      delete dsItr->second;
    dataStructMap.clear();
  }
  
  void 
  SolverAssembly::readInput(Lucee::LuaTable& tbl)
  {
// get stream for logging
    Lucee::Logger& logger = Lucee::Logger::get("lucee.console");
    Lucee::LogStream infoStrm = logger.getInfoStream();

// create all grids in assembly
    std::vector<std::string> gridNms = tbl.getNamesOfType("Grid");
    for (unsigned i=0; i<gridNms.size(); ++i)
    {
      std::string gnm = gridNms[i];
      Lucee::LuaTable gtbl = tbl.getTable(gnm);
      std::string kind = gtbl.getKind();
      infoStrm << "Setting up grid '" << gnm << "' of kind '"
               << kind << "'" << std::endl;
      std::auto_ptr<Lucee::GenericFactory<Lucee::GridIfc> > gfact(
        Lucee::ObjCreator<Lucee::GenericFactory<Lucee::GridIfc> >::getNew(kind));
      gfact->readInput(gtbl); // read input from grid block
      gridMap[gnm] = gfact->create(*this); // create grid
    }

// create all data-structs in assembly
    std::vector<std::string> dsNms = tbl.getNamesOfType("DataStruct");
    for (unsigned i=0; i<dsNms.size(); ++i)
    {
      std::string dsnm = dsNms[i];
      Lucee::LuaTable dstbl = tbl.getTable(dsnm);
      std::string kind = dstbl.getKind();
      infoStrm << "Setting up data-structure '" << dsnm << "' of kind '"
               << kind << "'" << std::endl;
      std::auto_ptr<Lucee::GenericFactory<Lucee::DataStructIfc> > dsfact(
        Lucee::ObjCreator<Lucee::GenericFactory<Lucee::DataStructIfc> >::getNew(kind));
      dsfact->readInput(dstbl); // read input from grid block
      dataStructMap[dsnm] = dsfact->create(*this); // create grid
    }
  }

  void
  SolverAssembly::buildData()
  {
  }

  void 
  SolverAssembly::buildAlgorithms()
  {
  }

  void
  SolverAssembly::initialize()
  {
  }
  
  void 
  SolverAssembly::writeToFile(const std::string& baseName, unsigned d)
  {
// create HDF5 for storing data
    Lucee::Hdf5Io io(0, 0);
    std::ostringstream fn;
    fn << baseName << "_" << d << ".h5";
    Lucee::IoNodeType fNode = io.createFile(fn.str());

// write out all grids
    std::map<std::string, Lucee::GridIfc*>::iterator gItr;
    for (gItr = gridMap.begin(); gItr != gridMap.end(); ++gItr)
      gItr->second->writeToFile(io, fNode, gItr->first);
  }

  int
  SolverAssembly::advance(double t)
  {
    return 0;
  }

  void
  SolverAssembly::restoreFromFile(const std::string& baseName)
  {
  }

  void
  SolverAssembly::finalize()
  {
  }
}
