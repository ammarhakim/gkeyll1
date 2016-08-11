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
 * Set direction to apply BC.
 *
 * @param d Direction in which BC is being applied.
 */
      void setDir(unsigned d) { dir = d; }

/**
 * Set edge to which BC is being applied.
 *
 * @param e Edge to which BC is being applied.
 */
      void setEdge(unsigned e) { edge = e; }

/** 
 * Apply boundary conditions. The variable 'qin' represents data in
 * the skin cell (first/last interior cell) while 'qbc' represents
 * data in the ghost cells. Also passed is the data in the cell "left"
 * of the skin cell, to allow for linear extrapolation BCs if
 * needed. Usually, this can be ignored.
 *
 * @param tm Time at which BC should be applied.
 * @param loc Coordinate at which BC is applied.
 * @param idx Grid index at which BC is applied.
 * @param c Coordinate system to use.
 * @param qin1 Data in cell "left" of skin cell.
 * @param qin Data in skin cell.
 * @param qbc Data in ghost cell.
 */
      virtual void applyBc(double tm, const double loc[3], const int *idx,
        const Lucee::RectCoordSys& c, 
        const Lucee::ConstFieldPtr<double>& qin1, const Lucee::ConstFieldPtr<double>& qin,
        Lucee::FieldPtr<double>& qbc) = 0;

    protected:
/**
 * Set components: this can be used by derived classes to set the
 * components array.
 *
 * @param comps Component list to set.
 */
      void setComponents(const std::vector<unsigned>& comps);

/**
 * Direction in which BC is being applied.
 *
 * @return Direction in which BC is being applied.
 */
      unsigned getDir() const { return dir; }

/**
 * Edge to which BC is being applied.
 *
 * @return Edge to which BC is being applied.
 */
      unsigned getEdge() const { return edge; }

    private:
/** Direction to apply boundary condtion */
      unsigned dir;
/** Edge to apply boundary condition */
      unsigned edge;
/** Components to apply to */
      std::vector<unsigned> components;
  };
}

#endif // LC_BOUNDARY_CONDITION_H
