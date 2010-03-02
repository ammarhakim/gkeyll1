/**
 * @file	LcSolverIfc.h
 *
 * @brief	Interface class for Lucee solvers.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_SOLVER_IFC_H
#define LC_SOLVER_IFC_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// std includes
#include <string>

namespace Lucee
{
/**
 * Base class for all solvers in Lucee. This class defines the basic
 * interface which solvers must support. Each solver has three
 * distinct types of interfaces (a) bootstrap, (b) time-step (c) I/O.
 *
 * The bootstrap interface methods allow solvers to initialize
 * themselves. The time-step interface methods expose methods to
 * advance solution in time. The I/O interface allow writing
 * simulation data to output files and restarting simulations from
 * previous output.
 */
  class SolverIfc
  {
    public:
/**
 * Bootstrap method: Read input from specified stream.
 */
      virtual void readInput() = 0;

/**
 * Bootstrap method: Allocate data for solver.
 */
      virtual void buildData() = 0;

/**
 * Initialize algorithms needed for solver.
 */
      virtual void buildAlgorithms() = 0;

/**
 * Initialize solver, i.e. setup initial conditions. At the end of
 * this call, the solver should be ready for evolving the solution.
 */
      virtual void initialize() = 0;

/**
 * Set the current time for solver. This is the time at which the
 * present state of the solver is valid.
 *
 * @param tm Current time.
 */
      virtual void setCurrTime(double tm) { currTime = tm; }

/**
 * Get the current time for solver. This is the time at which the
 * present state of the solver is valid.
 *
 * @retrun Current time.
 */
      virtual double getCurrTime() const { return currTime; }

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
 */
      virtual void writeToFile(const std::string& baseName) const = 0;

/**
 * Restore solver data from file. This is called instead of the
 * initialize() method if the simulation is being restarted.
 *
 * @param baseName Base name of input files. This should serves as a
 *   prefix for all input files.
 */
      virtual void restoreFromFile(const std::string& baseName) const = 0;

/**
 * Finalize solver: free resources, deallocate memory, close files
 * etc.
 */
      virtual void finalize() = 0;
      
    private:
/** Solver name */
      std::string name;
/** Current time */
      double currTime;
  };
}

#endif // LC_SOLVER_IFC_H
