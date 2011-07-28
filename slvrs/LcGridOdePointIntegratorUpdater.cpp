/**
 * @file	LcGridOdePointIntegratorUpdater.cpp
 *
 * @brief	Updater to integrate ODEs on a grid.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGridOdePointIntegratorUpdater.h>

namespace Lucee
{
  template <> const char *GridOdePointIntegratorUpdater<1>::id = "GridOdePointIntegrator1D";
  template <> const char *GridOdePointIntegratorUpdater<2>::id = "GridOdePointIntegrator2D";
  template <> const char *GridOdePointIntegratorUpdater<3>::id = "GridOdePointIntegrator3D";

  template <unsigned NDIM>
  GridOdePointIntegratorUpdater<NDIM>::GridOdePointIntegratorUpdater()
    : Lucee::UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  GridOdePointIntegratorUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method first
    Lucee::UpdaterIfc::readInput(tbl);
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  GridOdePointIntegratorUpdater<NDIM>::update(double t)
  {
    double dt = t-this->getCurrTime();

    return Lucee::UpdaterStatus(true, dt);
  }

  template <unsigned NDIM>
  void
  GridOdePointIntegratorUpdater<NDIM>::declareTypes()
  {
  }

// instantiations
  template class GridOdePointIntegratorUpdater<1>;
  template class GridOdePointIntegratorUpdater<2>;
  template class GridOdePointIntegratorUpdater<3>;
}
