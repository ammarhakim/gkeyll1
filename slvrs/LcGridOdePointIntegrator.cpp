/**
 * @file   LcGridOdePointIntegrator.cpp
 *
 * @brief   ODE integrator over complete grid.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGridOdePointIntegrator.h>

namespace Lucee
{
  template <unsigned NDIM>
  GridOdePointIntegrator<NDIM>::GridOdePointIntegrator(const Lucee::StructuredGridBase<NDIM>& grid)
    : GridOdeIntegrator<NDIM>(grid, 1)
  {
  }

  template <unsigned NDIM>
  void
  GridOdePointIntegrator<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
  }

  template <unsigned NDIM>
  void
  GridOdePointIntegrator<NDIM>::integrate(double t, Lucee::Field<NDIM, double>& sol)
  {
// fetch input field
    const Lucee::Field<NDIM, double>& inp = this->getIn(0);
// time-step to use
    double dt = t-this->getCurrTime();
// update using RK4 scheme
    rk4(dt, inp, sol);
  }

  template <unsigned NDIM>
  void
  GridOdePointIntegrator<NDIM>::rk4(double dt, const Lucee::Field<NDIM, double>& inp,
    Lucee::Field<NDIM, double>& sol)
  {
  }

// instantiations
  template class GridOdePointIntegrator<1>;
  template class GridOdePointIntegrator<2>;
  template class GridOdePointIntegrator<3>;
}
