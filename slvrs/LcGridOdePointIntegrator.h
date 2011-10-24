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
#include <LcPointSourceIfc.h>

// std include
#include <vector>

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
 * Delete object.
 */
      virtual ~GridOdePointIntegrator();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Integrate ODEs to time "t".
 *
 * @param t0 Initial time
 * @param t1 Final time at which solution is required.
 * @param sol On input, initial solution at t0. On output, solution at t1.
 */
      void integrate(double t0, double t1, Lucee::Field<NDIM, double>& sol);

    private:
/** Integration scheme to use */
      unsigned odeScheme;
/** Point sources to use as RHS in integrator */
      std::vector<Lucee::PointSourceIfc*> rhs;
/** Vector for use in computing sources */
      std::vector<double> ts;
/** Switch to indicate integration method */
      unsigned intMethod;

/**
 * Integrate ODEs to time "t" using Runge-Kutta 4th order method.
 *
 * @param t0 Initial time
 * @param dt Time-step to use
 * @param sol On output, solution of ODEs.
 */
      void rk4(double t0, double dt, Lucee::Field<NDIM, double>& sol);

/**
 * Integrate ODEs to time "t" using semi-implicit method.
 *
 * @param t0 Initial time
 * @param dt Time-step to use
 * @param sol On output, solution of ODEs.
 */
      void semiImplicit(double t0, double dt, Lucee::Field<NDIM, double>& sol);

/**
 * Compute sources, summing up contributions from each RHS term.
 *
 * @param tm Time at which source is requested.
 * @param xc Coordinates at which source is needed
 * @param inp Inputs for which source is needed
 * @param src On output, sources.
 */
      void calcSource(double tm, const double xc[3], const double *inp, std::vector<double>& src);

/**
 * Compute source jacobian, summing up contributions from each RHS term.
 *
 * @param tm Time at which source is requested.
 * @param xc Coordinates at which source is needed
 * @param inp Inputs for which source is needed
 * @param srcJac On output, source Jacobian
 */
      void calcSourceJac(double tm, const double xc[3], const double *inp, 
        Lucee::Matrix<double>& srcJac);
  };
}

#endif // LC_GRID_ODE_POINT_INTEGRATOR_H
