/**
 * @file	LcBoundaryCondition.h
 *
 * @brief	Base class for applying boundary condition.
 */

#ifndef LC_BOUNDARY_CONDITION_H
#define LC_BOUNDARY_CONDITION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcConstFieldPtr.h>
#include <LcFieldPtr.h>
#include <LcRectCoordSys.h>

namespace Lucee
{
/** Interface class for applying boundary conditions */
  class BoundaryCondition : public Lucee::BasicObj
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Default ctor */
      BoundaryCondition();

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
 * @param qout Data in ghost cell.
 */
      virtual void applyBc(const Lucee::RectCoordSys& c, const Lucee::ConstFieldPtr<double>& qin, Lucee::FieldPtr<double>& qbc) = 0;
  };
}

#endif // LC_BOUNDARY_CONDITION_H
