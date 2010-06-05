/**
 * @file	LcSolverAssembly.h
 *
 * @brief	Solver that assembles updaters to create a simulation.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_SOLVER_ASSEMBLY_H
#define LC_SOLVER_ASSEMBLY_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGridBase.h>
#include <LcSolverIfc.h>

// std includes
#include <map>
#include <string>

namespace Lucee
{
/**
 * A general purpose solver that assembles a particular type of
 * solvers (called Updaters) to perform a simulation. Grids and
 * data-structures, in addition to updaters, can be used in the
 * simulation. Various time-stepping schemes can be used to advance
 * the solution in time.
 */
  class SolverAssembly : public Lucee::SolverIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create a new assembly object 
 */
      SolverAssembly();

/**
 * Destroy the object
 */
      virtual ~SolverAssembly();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Bootstrap method: Allocate data for solver.
 */
      virtual void buildData();

/**
 * Initialize algorithms needed for solver.
 */
      virtual void buildAlgorithms();

/**
 * Initialize solver, i.e. setup initial conditions. At the end of
 * this call, the solver should be ready for evolving the solution.
 */
      virtual void initialize();

/**
 * Advance the solution to specified time. Solvers that do not have a
 * concept of time should ignore the time parameter.
 *
 * @param t Time to advance the solution to.
 * @return Status of solution.
 */
      virtual int advance(double t);

/**
 * Write solver data to file.
 *
 * @param baseName Base name of output files. This should serve as a
 *   prefix for all output files.
 * @param d Dump number.
 */
      virtual void writeToFile(const std::string& baseName, unsigned d);

/**
 * Restore solver data from file. This is called instead of the
 * initialize() method if the simulation is being restarted.
 *
 * @param baseName Base name of input files. This should serves as a
 *   prefix for all input files.
 */
      virtual void restoreFromFile(const std::string& baseName);

/**
 * Finalize solver: free resources, deallocate memory, close files
 * etc.
 */
      virtual void finalize();

    private:
/** Map of grids */
      std::map<std::string, Lucee::GridBase*> gridMap;
/** Map of fields */

/** Map of field-aliases */

/** Map of updaters */

  };
}

#endif // LC_SOLVER_ASSEMBLY_H
