/**
 * @file	LcDistFuncReflectionBcUpdater.h
 *
 * @brief	Applies particle refection BCs to distribution function
 */

#ifndef LC_DIST_FUNC_REFLECTION_BS_UPDATER_H
#define LC_DIST_FUNC_REFLECTION_BS_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcMatrix.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>
#include <LcVector.h>

namespace Lucee
{
/**
 * Applies particle refection BCs to distribution function
 */
  class DistFuncReflectionBcUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new projection updater */
      DistFuncReflectionBcUpdater();

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
 * Advance the solution to specified time. Updaters that do not have a
 * concept of time should ignore the time parameter.
 *
 * @param t Time to advance the solution to.
 * @return Status of updater.
 */
      Lucee::UpdaterStatus update(double t);

/**
 * Declare the types of input and output variables accepted by this
 * updater. This must be provided by the derived classes for the
 * type-checking to pass. Inside the implementation of this method the
 * derived class must make a sequence of appendInpVarType() and
 * appendOutVarType() calls to declare the input/output data structure
 * types.
 */
      void declareTypes();

/**
 * Method that performs registration of Lua functions.
 *
 * @param lfm Lua function map object.
 */
      static void appendLuaCallableMethods(Lucee::LuaFuncMap& lfm);

    private:
/** Pointer to nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<2> *nodalBasis;
/** Flag to indicate if BCs should be applied to left edge */
      bool applyLeftEdge;
/** Flag to indicate if BCs should be applied to right edge */
      bool applyRightEdge;
/** Mapping for 180 degree rotations for upper edge */
      std::vector<unsigned> rotMapRight;
/** Mapping for 180 degree rotations for lower edge */
      std::vector<unsigned> rotMapLeft;
/** Cutoff velocity */
      double cutOffVel;

/**
 * Lua callable method to set cut-off velocity.
 *
 * @param L Lua state to use.
 * @return number of output parameters.
 */
      static int luaSetCutOffVelocity(lua_State *L);

/**
 * Set cut-off velocity.
 *
 * @param cv Cut-off velocity.
 */
      void setCutOffVelocity(double cv)
      { 
        cutOffVel = cv;
      }
  };
}

#endif // LC_DIST_FUNC_REFLECTION_BS_UPDATER_H
