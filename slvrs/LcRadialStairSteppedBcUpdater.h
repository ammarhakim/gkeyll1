/**
 * @file	LcRadialStairSteppedBcUpdater.h
 *
 * @brief	Base class for boundary conditions for radial stair-stepped boundaries.
 */

#ifndef LC_RADIAL_STAIR_STEPPED_BC_UPDATER_H
#define LC_RADIAL_STAIR_STEPPED_BC_UPDATER_H

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
 * Base class for boundary conditions for radial stair-stepped
 * boundaries.
 */
  template <unsigned NDIM>
  class RadialStairSteppedBcUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create new updater.
 */
      RadialStairSteppedBcUpdater();

/**
 * Delete updater.
 */
      virtual ~RadialStairSteppedBcUpdater();

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
 * Write grid to specified HDF5 file.
 *
 * @param nm Name of file to write.
 */
      void write(const std::string& nm);

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
 * Lua callable method to set direction to update.
 *
 * @param L Lua state to use.
 * @return number of output parameters.
 */
      static int luaSetDir(lua_State *L);

/**
 * Lua callable method for writing out grid data to HDF5 file.
 *
 * @param L Lua state to use.
 * @return number of output parameters.
 */
      static int luaWrite(lua_State *L);


    private:
/** Direction to apply boundary conditon */
      unsigned bcDir;
/** Origin of stair stepped boundary */
      Lucee::Vec3<double> origin;
/** Radius of stair stepped boundary */
      double radius;
/** Components to be set constant */
      std::vector<unsigned> constComponents;
      std::vector<double> constValues;
/** Components to be radially copied */
      std::vector<unsigned> copyComponents;
/** factor to multiply while copying */
      std::vector<double> copyFact;
/** Components to be radially reflected */
      std::vector<unsigned> reflectComponents;
/** Components to be radially absorbed */
      std::vector<unsigned> absorbComponents;
/** Components to be transversely reflected */
      std::vector<unsigned> reflectTransComponents;
/** Components with zero radial components */
      std::vector<unsigned> zeroRadialComponents;
/** Components with zero transverse components */
      std::vector<unsigned> zeroTransComponents;
/** Components to be computed using a function */
      std::vector<unsigned> fieldFunctionComponents;
/** Input components for the field function BC*/
      std::vector<unsigned> fieldFunctionInpComponents;
/** Reference to Lua function for field function BC */
      int fnRef;
/** Flag if we pre-process single cell contribution before accumulating **/
      bool preProcess;
/** Reference to Lua function for pre-processing field function*/
      int preProcessFnRef;
/** Pointer to in/out field */
      Lucee::Field<NDIM, double> *inOut;

/**
 * Set direction to update.
 *
 * @param dir Direction to update
 */
      void setDir(unsigned dir) { bcDir = dir; }
  };
}

#endif // LC_STAIR_STEPPED_BC_UPDATER_H
