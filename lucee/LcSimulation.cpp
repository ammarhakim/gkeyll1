/**
 * @file	LcSimulation.cpp
 *
 * @brief	Top-level simulation class.
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
#include <LcObjCreator.h>
#include <LcSimulation.h>


namespace Lucee
{
  Simulation::Simulation()
    : Lucee::SolverIfc("Simulation")
  {
  }

  void
  Simulation::readInput(Lucee::LuaTable& tbl)
  {
    slvr = Lucee::ObjCreator<Lucee::SolverIfc>::getNew(tbl.getKind());
    slvr->readInput(tbl);
  }

  void
  Simulation::buildData()
  {
    slvr->buildData();
  }

  void
  Simulation::buildAlgorithms()
  {
    slvr->buildAlgorithms();
  }

  void
  Simulation::initialize()
  {
    slvr->initialize();
  }

  int
  Simulation::advance(double t)
  {

    return 0;
  }

  void
  Simulation::writeToFile(const std::string& baseName) const
  {
  }

  void
  Simulation::restoreFromFile(const std::string& baseName)
  {
  }

  void
  Simulation::finalize()
  {
  }
}
