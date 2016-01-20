/**
 * @file	LcUpdaterIfc.h
 *
 * @brief	Base class for updaters in Lucee.
 */

#ifndef LC_UPDATER_IFC_H
#define LC_UPDATER_IFC_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBasicObj.h>
#include <LcDataStructIfc.h>
#include <LcGridIfc.h>
#include <LcUpdaterStatus.h>

// std includes
#include <string>
#include <typeinfo>
#include <vector>

namespace Lucee
{
/**
 * An updater is analogous to a subroutine: it takes a set of
 * input/output DataStructIfc objects and advances them to a specified
 * time.
 *
 * The updater provides two sets of methods: initialization and
 * advancement. The initialization sequence is: 
 *
 * derived class construction methods if the updator is used directly
 * from code or readInput() if the updater is being used from input
 * file.
 *
 * declareTypes() to declare the input/output types the updater accepts.
 * 
 * setGrid() to set the grid on which the updater should run.
 *
 * To advancement sequence is:
 *
 * calls to setInpVars(), setOutVars() to set the input/output data
 * structures used by the updater.
 *
 * setCurrTime(tm) to set the current time (i.e. the time at which the
 * solution is already computed).
 *
 * update(t) to advance the solution to time t. The time-step to use
 * is dt = t-getCurrTime(). This method must return an instance of
 * UpdaterStatus class, indicating the status of the update() method.
 */
  class UpdaterIfc : public Lucee::BasicObj
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
 * Initialize solver, i.e. setup initial conditions. At the end of
 * this call, the solver should be ready for evolving the solution.
 */
      virtual void initialize();

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
      double getCurrTime() const;

/**
 * Advance the solution to specified time. Updaters that do not have a
 * concept of time should ignore the time parameter.
 *
 * @param t Time to advance the solution to.
 * @return Status of updater.
 */
      virtual Lucee::UpdaterStatus update(double t) = 0;

/**
 * Declare the types of input and output variables accepted by this
 * updater. This must be provided by the derived classes for the
 * type-checking to pass. Inside the implementation of this method the
 * derived class must make a sequence of appendInpVarType() and
 * appendOutVarType() calls to declare the input/output data structure
 * types.
 */
      virtual void declareTypes() = 0;

/**
 * Set the grid on which to run this updater.
 *
 * @param grd Grid to run this updater on.
 */
      void setGrid(const Lucee::GridIfc& grd);

/**
 * Set input data-structures.
 *
 * @param dsl List of input data structures.
 */
      void setInpVars(const std::vector<const Lucee::DataStructIfc*>& dsl);

/**
 * Set output data-structures.
 *
 * @param dsl List of output data structures.
 */
      void setOutVars(const std::vector<Lucee::DataStructIfc*>& dsl);

/**
 * Get number of input data structures passed to updater.
 *
 * @return number of input data structures.
 */
      unsigned getNumInpVars() const
      { 
        return inpVars.size();
      }

/**
 * Get number of output data structures passed to updater.
 *
 * @return number of output data structures.
 */
      unsigned getNumOutVars() const
      {
        return outVars.size();
      }

/**
 * Method that performs registration of Lua functions.
 *
 * @param lfm Lua function map object.
 */
      static void appendLuaCallableMethods(Lucee::LuaFuncMap& lfm);

/**
 * Lua callable method to set current time.
 *
 * @param L Lua state to use.
 * @return number of output parameters.
 */
      static int luaSetCurrTime(lua_State *L);

/**
 * Lua callable method to advance solver to given time.
 *
 * @param L Lua state to use.
 * @return number of output parameters.
 */
      static int luaAdvance(lua_State *L);

/**
 * Lua callable method to set input datastructures.
 *
 * @param L Lua state to use.
 * @return number of output parameters.
 */
      static int luaSetInpVars(lua_State *L);

/**
 * Lua callable method to set output datastructures.
 *
 * @param L Lua state to use.
 * @return number of output parameters.
 */
      static int luaSetOutVars(lua_State *L);

/**
 * Lua callable method to get wall-clock time for updater
 *
 * @param L Lua state to use.
 * @return number of output parameters.
 */
      static int luaGetTime(lua_State *L);

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
        if (grid == 0)
          throw Lucee::Except("UpdaterIfc::getGrid: grid pointer is not valid");
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
/** Current time */
      double currTime;
/** Grid on which updater should be applied */
      const Lucee::GridIfc *grid;
/** List of input data structures */
      std::vector<const Lucee::DataStructIfc*> inpVars;
/** List of output data structures */
      std::vector<Lucee::DataStructIfc*> outVars;
/** Cumulative wall-clock time spent in the advance method */
      double totAdvanceWallTime;
/** Cumulative CPU time spent in the advance method */
      double totAdvanceCpuTime;
  };
}

#endif // LC_UPDATER_IFC_H
