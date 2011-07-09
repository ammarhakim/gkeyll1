/**
 * @file	LcZeroTangentBoundaryCondition.h
 *
 * @brief	Class for applying zero-tangent BCs.
 */

#ifndef LC_ZERO_TANGENT_BOUNDARY_CONDITION_H
#define LC_ZERO_TANGENT_BOUNDARY_CONDITION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBoundaryCondition.h>

namespace Lucee
{
/** Class to apply a zero-tangent boundary condition to vectors */
  class ZeroTangentBoundaryCondition : public Lucee::BoundaryCondition
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
 * @param c Coordinate system to use.
 * @param qin Data in skin cell.
 * @param qbc Data in ghost cell.
 */
      void applyBc(const Lucee::RectCoordSys& c, const Lucee::ConstFieldPtr<double>& qin, Lucee::FieldPtr<double>& qbc);
  };
}

#endif // LC_ZERO_TANGENT_BOUNDARY_CONDITION_H
