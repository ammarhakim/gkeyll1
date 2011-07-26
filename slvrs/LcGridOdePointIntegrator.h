/**
 * @file   LcGridOdePointIntegrator.h
 *
 * @brief   ODE integrator over complete grid.
 */

#ifndef LC_GRID_ODE_POINT_INTEGRATOR_H
#define LC_GRID_ODE_POINT_INTEGRATOR_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGridOdeIntegrator.h>

namespace Lucee
{
/**
 * Class for ODE integrator for functions defined on each cell of a
 * grid.
 */
  template <unsigned NDIM>
  class GridOdePointIntegrator : public GridOdeIntegrator<NDIM>
  {
    public:
/**
 * Create ODE integrator to integrate equations on specified grid.
 *
 * @param grid Grid to integrate ODEs on.
 */
      GridOdePointIntegrator(const Lucee::StructuredGridBase<NDIM>& grid);

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Integrate ODEs to time "t".
 *
 * @param t Time to integrate ODEs to.
 * @param sol On output, solution of ODEs.
 */
      void integrate(double t, Lucee::Field<NDIM, double>& sol);

    private:
/** Integration scheme to use */
      unsigned odeScheme;

/**
 * Integrate ODEs to time "t".
 *
 * @param t Time to integrate ODEs to.
 * @param inp Input field.
 * @param sol On output, solution of ODEs.
 */
      virtual void rk4(double dt, const Lucee::Field<NDIM, double>& inp,
        Lucee::Field<NDIM, double>& sol);
  };
}

#endif // LC_GRID_ODE_POINT_INTEGRATOR_H
