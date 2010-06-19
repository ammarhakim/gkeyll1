/**
 * @file	LcUpdaterIfc.h
 *
 * @brief	Base class for updaters in Lucee.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_UPDATER_IFC_H
#define LC_UPDATER_IFC_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDataStructIfc.h>
#include <LcGridIfc.h>
#include <LcHdf5Io.h>
#include <LcSolverIfc.h>

// std includes
#include <map>
#include <string>
#include <vector>

// forward declare SolverAssembly class
class Lucee::SolverAssembly;

namespace Lucee
{
/**
 * An updater is analogous to a subroutine: it takes a set of
 * input/output DataStructIfc objects and advances them to a specified
 * time.
 */
  class UpdaterIfc : public Lucee::SolverIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create a new updater object 
 */
      UpdaterIfc();

/**
 * Destroy the object
 */
      virtual ~UpdaterIfc();

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
 * Advance the solution to specified time. Updaters that do not have a
 * concept of time should ignore the time parameter.
 *
 * @param t Time to advance the solution to.
 * @return Status of solution.
 */
      virtual int advance(double t) = 0;

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

    protected:
/**
 * Get input dataStruct at specified location.
 *
 * @param loc Location in input dataStruct list.
 * @return const reference to dataStruct.
 */
      template <typename DS>
      const DS& getInp(unsigned loc) const
      {
        return parent->template
          getConstDataStruct<DS>(inpVarNames[loc]);
      }

/**
 * Get output dataStruct at specified location.
 *
 * @param loc Location in output dataStruct list.
 * @return reference to dataStruct.
 */
      template <typename DS>
      DS& getOut(unsigned loc)
      {
        return parent->template
          getDataStruct<DS>(inpVarNames[loc]);
      }

    private:
/** List of input data structure names */
      std::vector<std::string> inpVarNames;
/** List of output data structure names */
      std::vector<std::string> outVarNames;
/** Name of domain on which updater should be applied */
      std::string onGrid;
/** Pointer to containing solver assembly */
      Lucee::SolverIfc *parent;
  };
}

#endif // LC_UPDATER_IFC_H
