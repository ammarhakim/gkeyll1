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
    std::map<std::string, Lucee::GridBase*>::iterator gItr;
    for (gItr = gridMap.begin(); gItr != gridMap.end(); ++gItr)
      delete gItr->second;
    gridMap.clear();
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
      std::auto_ptr<Lucee::GenericFactory<Lucee::GridBase> > gfact(
        Lucee::ObjCreator<Lucee::GenericFactory<Lucee::GridBase> >::getNew(kind));
      gfact->readInput(gtbl); // read input from grid block
      gridMap[gnm] = gfact->create(); // create grid
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
    std::map<std::string, Lucee::GridBase*>::iterator gItr;
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
