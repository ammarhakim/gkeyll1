/**
 * @file	LcSimulation.h
 *
 * @brief	Top-level simulation class.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_SIMULATION_H
#define LC_SIMULATION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCmdLineArgs.h>
#include <LcSolverIfc.h>

// std includes
#include <string>

namespace Lucee
{
/**
 * Top-level class to drive all simulations.
 */
  class Simulation : public Lucee::SolverIfc
  {
    public:
/**
 * Create a new simulation object.
 */
      Simulation();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      void readInput(Lucee::LuaTable& tbl);

/**
 * Bootstrap method: Allocate data for solver.
 */
      void buildData();

/**
 * Initialize algorithms needed for solver.
 */
      void buildAlgorithms();

/**
 * Initialize solver, i.e. setup initial conditions. At the end of
 * this call, the solver should be ready for evolving the solution.
 */
      void initialize();

/**
 * Advance the solution to specified time. Solvers that do not have a
 * concept of time should ignore the time parameter.
 *
 * @param t Time to advance the solution to.
 * @return Status of solution.
 */
      int advance(double t);

/**
 * Write solver data to file.
 *
 * @param baseName Base name of output files. This should serve as a
 *   prefix for all output files.
 */
      void writeToFile(const std::string& baseName) const;

/**
 * Restore solver data from file. This is called instead of the
 * initialize() method if the simulation is being restarted.
 *
 * @param baseName Base name of input files. This should serves as a
 *   prefix for all input files.
 */
      void restoreFromFile(const std::string& baseName);

/**
 * Finalize solver: free resources, deallocate memory, close files
 * etc.
 */
      void finalize();

    private:
/** Top-level solver */
      Lucee::SolverIfc *slvr;
  };
}

#endif // LC_SIMULATION_H
