/**
 * @file	LcSolverIfc.h
 *
 * @brief	Interface class for Lucee solvers.
 */

#ifndef LC_SOLVER_IFC_H
#define LC_SOLVER_IFC_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLuaTable.h>
#include <LcBasicObj.h>

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
 * previous output. The order in the bootstrap methods is called
 * readInput().
 *
 * Next, if the simulation is a restart, the method restoreFromFile()
 * is called. Otherwise, if the simulation is not a restart, the
 * method initialize() is called.
 *
 * The advance phase methods are called in the order setCurrTime(told)
 * and then advance(t). The solver should advance itself by time
 * t-told. The return code of the advance() method tells Lucee if
 * solver failed or passed and what should be done next.
 *
 * The writeToFile() method is called for the solver to write out its
 * data to a file.
 *
 * The finalize() method is called to shut-down the solver.
 */
  class SolverIfc : public Lucee::BasicObj
  {
    public:
/** Class id: this is used by the registration system */
      static const char *id;

/**
 * Create a new solver object with given name.
 *
 * @param nm Name of solver.
 */
      SolverIfc(const std::string& nm);

/**
 * Destroy solver object.
 */
      virtual ~SolverIfc();

/**
 * Set number of output frames to write. This does not include the
 * first and last frames.
 */
      void setNumOutFrames(unsigned n);

/**
 * Get number of output frames to write. This does not include the
 * first and last frames.
 *
 * @return number of output frames.
 */
      unsigned getNumOutFrames(unsigned n) const { return numOutFrames; }

/**
 * Set the current time for solver. This is the time at which the
 * present state of the solver is valid.
 *
 * @param tm Current time.
 */
      void setCurrTime(double tm);

/**
 * Get the current time for solver. This is the time at which the
 * present state of the solver is valid.
 *
 * @return Current time.
 */
      double getCurrTime() const { return currTime; }

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Initialize solver, i.e. setup initial conditions. At the end of
 * this call, the solver should be ready for evolving the solution.
 */
      virtual void initialize() = 0;

/**
 * Advance the solution to specified time. Solvers that do not have a
 * concept of time should ignore the time parameter. The solver should
 * return 0 if it fails and 1 it passes.
 *
 * @param t Time to advance the solution to.
 * @return Status of solution: 1 for pass, 0 for fail.
 */
      virtual int advance(double t) = 0;

/**
 * Write solver data to file.
 *
 * @param baseName Base name of output files. This should serve as a
 *   prefix for all output files.
 * @param d Dump number.
 */
      virtual void writeToFile(const std::string& baseName, unsigned d) = 0;

/**
 * Restore solver data from file. This is called instead of the
 * initialize() method if the simulation is being restarted.
 *
 * @param baseName Base name of input files. This should serves as a
 *   prefix for all input files.
 */
      virtual void restoreFromFile(const std::string& baseName) = 0;

/**
 * Finalize solver: free resources, deallocate memory, close files
 * etc.
 */
      virtual void finalize() = 0;

/**
 * Method that performs registration of Lua functions.
 *
 * @param lfm Lua function map object.
 */
      static void appendLuaCallableMethods(Lucee::LuaFuncMap& lfm);

/**
 * Lua callable method to advance solver to given time.
 *
 * @param L Lua state to use.
 * @return number of output parameters.
 */
      static int luaAdvance(lua_State *L);

/**
 * Lua callable method to initialize solver to given time.
 *
 * @param L Lua state to use.
 * @return number of output parameters.
 */
      static int luaInitialize(lua_State *L);

/**
 * Lua callable method to write data to file.
 *
 * @param L Lua state to use.
 * @return number of output parameters.
 */
      static int luaWrite(lua_State *L);

    private:
/** Current time */
      double currTime;
/** Number of output frames to write */
      unsigned numOutFrames;
  };
}

#endif // LC_SOLVER_IFC_H
