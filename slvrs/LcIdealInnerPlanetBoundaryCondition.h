/**
 * @file	LcIdealInnerPlanetBoundaryCondition.h
 *
 * @brief	Class for applying BCs to inner planet surface based on ideal-MHD
 */

#ifndef LC_IDEAL_INNER_PLANET_BOUNDARY_CONDITION_H
#define LC_IDEAL_INNER_PLANET_BOUNDARY_CONDITION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBoundaryCondition.h>

namespace Lucee
{
/** Class to apply BCs to inner planet surface based on ideal-MHD */
  class IdealInnerPlanetBoundaryCondition : public Lucee::BoundaryCondition
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
 * @param idx Grid index at which source is requested.
 * @param c Coordinate system to use.
 * @param qin Data in skin cell.
 * @param qbc Data in ghost cell.
 */
      void applyBc(double tm, const double loc[3], const int *idx,
        const Lucee::RectCoordSys& c, const Lucee::ConstFieldPtr<double>& qin1,
        const Lucee::ConstFieldPtr<double>& qin, Lucee::FieldPtr<double>& qbc);

    private:
  };
}

#endif // LC_IDEAL_INNER_PLANET_BOUNDARY_CONDITION_H
