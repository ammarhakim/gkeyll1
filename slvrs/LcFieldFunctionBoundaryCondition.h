/**
 * @file	LcFieldFunctionBoundaryCondition.h
 *
 * @brief	Class for applying boundary conditions specified via a Lua function.
 */

#ifndef LC_FIELD_FUNCTION_BOUNDARY_CONDITION_H
#define LC_FIELD_FUNCTION_BOUNDARY_CONDITION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBoundaryCondition.h>

namespace Lucee
{
/** Interface class for applying boundary conditions */
  class FieldFunctionBoundaryCondition : public Lucee::BoundaryCondition
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/** 
 * Apply boundary conditions. The variable 'qin' represents data in
 * the skin cell (first/last interior cell) while 'qbc' represents
 * data in the ghost cells.
 *
 * @param tm Time at which source is requested.
 * @param loc Coordinate at which source is requested.
 * @param c Coordinate system to use.
 * @param qin Data in skin cell.
 * @param qbc Data in ghost cell.
 */
      void applyBc(double tm, const double loc[3], const Lucee::RectCoordSys& c, 
        const Lucee::ConstFieldPtr<double>& qin, Lucee::FieldPtr<double>& qbc);

    private:
/** Reference to Lua function */
      int fnRef;
/** Input components to send to Lua function */
      std::vector<unsigned> inpComponents;
  };
}

#endif // LC_FIELD_FUNCTION_BOUNDARY_CONDITION_H
