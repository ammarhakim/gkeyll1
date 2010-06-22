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
#include <string>
#include <typeinfo>
#include <vector>

namespace Lucee
{
// forward declare SolverAssembly class
  class SolverAssembly;

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
 * @return Status of updater: 0 for pass, 1 otherwise.
 */
      virtual int advance(double t) = 0;

/**
 * Get a suggested time-step. This method is called after advance()
 * method is called. This method is called even if advance() fails
 * (i.e. return 0). If advance() returns 0 this step may be retaken
 * depending on the time-stepping mode used. If advance() returns 1
 * the next step may be larger than one supplied.
 *
 * @param suggested time-step.
 */
      virtual double getSuggestedDt() = 0;

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

/**
 * Set list of input variables for use in updater.
 *
 * @param nms list on input variables.
 */
      void setInpVarNames(const std::vector<std::string>& nms);

/**
 * Set list of output variables for use in updater.
 *
 * @param nms list on output variables.
 */
      void setOutVarNames(const std::vector<std::string>& nms);

/**
 * Get actual number of input data structures passed to updater.
 *
 * @return actual number of input data structures.
 */
      unsigned getActualNumInpVars() const
      { return inpVarNames.size(); }

/**
 * Get actual number of output data structures passed to updater.
 *
 * @return actual number of output data structures.
 */
      unsigned getActualNumOutVars() const
      { return outVarNames.size(); }

/**
 * Type-check if the supplied list of input/output data structures are
 * of the proper type.
 *
 * @param inp Names of input data structures.
 * @param out Names of output data structures.
 * 
 * @return true if type-checking passes, fail otherwise.
 */
      bool typeCheck(const std::vector<std::string>& inp,
        const std::vector<std::string>& out) const;

    protected:
/**
 * Set type information for an input data structure. This method must
 * be called for each input variable in the order they are expected.
 *
 * @param type Typeid of input variable.
 */
      void appendInpVarType(const std::type_info& type);

/**
 * Set type information for an output data structure. This method must
 * be called for each output variable in the order they are expected.
 *
 * @param type Typeid of output variable.
 */
      void appendOutVarType(const std::type_info& type);

/**
 * Set type information for the last input data structure. This method
 * should only be called if there are a variable number of input data
 * structures expected by the updater.
 *
 * @param type Typeid of last input variable.
 */
      void setLastInpVarType(const std::type_info& type);

/**
 * Set type information for the last output data structure. This
 * method should only be called if there are a variable number of
 * output data structures expected by the updater.
 *
 * @param type Typeid of last output variable.
 */
      void setLastOutVarType(const std::type_info& type);

/**
 * Get grid on which updater should be applied.
 *
 * @return reference to grid.
 */
      template <typename G>
      const G&
      getGrid() const
      {
        return parent->template getConstGrid<G>(onGrid);
      }

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
          getDataStruct<DS>(outVarNames[loc]);
      }

    private:
/**
 * Input/output variable typeids.
 */
      struct VarTypeIds
      {
/** Map of variable locations to typeids */
          std::vector<const std::type_info*> varTypes;
/** Last variable type (all others after this are assumed to be this) */
          const std::type_info* lastVarType;
      };

/** List of input data structure names */
      std::vector<std::string> inpVarNames;
/** List of output data structure names */
      std::vector<std::string> outVarNames;
/** Name of domain on which updater should be applied */
      std::string onGrid;
/** Pointer to containing solver assembly */
      Lucee::SolverIfc *parent;
/** Input variable types */
      VarTypeIds inpVarTypes;
/** Output variable types */
      VarTypeIds outVarTypes;
  };
}

#endif // LC_UPDATER_IFC_H
