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

// std includes
#include <vector>

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
 * Return component given index.
 *
 * @param k Component index.
 * @return Corresponding component.
 */
      unsigned component(unsigned k) const 
      { return components[k]; }

/**
 * Number of components.
 *
 * @return Number of componenents.
 */
      unsigned numComponents() const
      { return components.size(); }

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
      virtual void applyBc(double tm, const double loc[3], const int *idx,
        const Lucee::RectCoordSys& c, 
        const Lucee::ConstFieldPtr<double>& qin, Lucee::FieldPtr<double>& qbc) = 0;

    protected:
/**
 * Set components: this can be used by derived classes to set the
 * components array.
 *
 * @param comps Component list to set.
 */
      void setComponents(const std::vector<unsigned>& comps);

    private:
/** Components to apply to */
      std::vector<unsigned> components;
  };
}

#endif // LC_BOUNDARY_CONDITION_H
