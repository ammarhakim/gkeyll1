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
 *
 * When updaters are used directly, they must be first constructed,
 * either using the constructor or using the set of provided
 * methods. and then the following methods *must* be called in the
 * specified order: declareTypes(), setGrid(), initialize(). Updaters
 * may not work if these methods are not called in this order. If the
 * setGrid() method is called then the initialize() method must be
 * called before using the updater.
 *
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
 * By default this method will return the largest possible
 * time-step. This may not always be desirable and derived classes
 * should override this method.
 *
 * @param suggested time-step.
 */
      virtual double getSuggestedDt();

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
 * Declare the types of input and output variables accepted by this
 * updater. This must be provided by the derived classes for the
 * type-checking to pass.
 */
      virtual void declareTypes() = 0;

/**
 * Set the grid on which to run this updater.
 *
 * @param grd Grid to run this updater on.
 */
      void setGrid(const Lucee::GridIfc& grd);

/**
 * Set input data-structure at specified location.
 *
 * @param loc Location in input data structure list.
 * @param ds Data structure.
 */
      void setInpVar(unsigned loc, const Lucee::DataStructIfc& ds);

/**
 * Set output data-structure at specified location.
 *
 * @param loc Location in output data structure list.
 * @param ds Data structure.
 */
      void setOutVar(unsigned loc, Lucee::DataStructIfc& ds);

/**
 * Get number of input data structures passed to updater.
 *
 * @return number of input data structures.
 */
      unsigned getNumInpVars() const
      { 
        return inpVarTypes.varTypes.size();
      }

/**
 * Get number of output data structures passed to updater.
 *
 * @return number of output data structures.
 */
      unsigned getNumOutVars() const
      {
        return outVarTypes.varTypes.size();
      }

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
        return dynamic_cast<const G&>(*grid);
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
        return dynamic_cast<const DS&>(*inpVars[loc]);
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
        return dynamic_cast<DS&>(*outVars[loc]);
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

/** Input variable types */
      VarTypeIds inpVarTypes;
/** Output variable types */
      VarTypeIds outVarTypes;

/** Grid on which updater should be applied */
      const Lucee::GridIfc *grid;
/** List of input data structures */
      std::vector<const Lucee::DataStructIfc*> inpVars;
/** List of output data structures */
      std::vector<Lucee::DataStructIfc*> outVars;
  };
}

#endif // LC_UPDATER_IFC_H
