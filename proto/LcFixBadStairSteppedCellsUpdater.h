/**
 * @file	LcFixBadStairSteppedCellsUpdater.h
 *
 * @brief	Detect and fix degenerate cells in a stair-stepped mesh.
 */

#ifndef LC_FIX_BAD_STAIR_STEPPED_CELLS_UPDATER_H
#define LC_FIX_BAD_STAIR_STEPPED_CELLS_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBoundaryCondition.h>
#include <LcUpdaterIfc.h>

namespace Lucee
{
/**
 * Class to detect and fix degenerate cells in a stair-stepped mesh.
 */
  template <unsigned NDIM>
  class FixBadStairSteppedCellsUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create new updater.
 */
      FixBadStairSteppedCellsUpdater();

/**
 * Delete updater.
 */
      virtual ~FixBadStairSteppedCellsUpdater();

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

/**
 * Lua callable method to get number of bad cells
 *
 * @param L Lua state to use.
 * @return number of output parameters.
 */
      static int luaGetNumBadCells(lua_State *L);

    private:
/** Number of bad cells */
      unsigned numBadCells;
  };
}

#endif // LC_FIX_BAD_STAIR_STEPPED_CELLS_UPDATER_H
