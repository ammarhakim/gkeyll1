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
#include <LcStructuredGridBase.h>

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
  GridOdePointIntegratorUpdater<NDIM>::~GridOdePointIntegratorUpdater()
  {
    delete integrator;
  }

  template <unsigned NDIM>
  void
  GridOdePointIntegratorUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method first
    Lucee::UpdaterIfc::readInput(tbl);
// create integrator object
    const Lucee::StructuredGridBase<NDIM>& grd 
      = this->template getGrid<Lucee::StructuredGridBase<NDIM> >();
    integrator = new Lucee::GridOdePointIntegrator<NDIM>(grd);
// setup integrator
    integrator->readInput(tbl);
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  GridOdePointIntegratorUpdater<NDIM>::update(double t)
  {
// fetch field to update
    Lucee::Field<NDIM, double>& fld = this->getOut<Lucee::Field<NDIM, double> >(0);
// integrate ODEs
    integrator->integrate(this->getCurrTime(), t, fld);

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  GridOdePointIntegratorUpdater<NDIM>::declareTypes()
  {
// expect one output field
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class GridOdePointIntegratorUpdater<1>;
  template class GridOdePointIntegratorUpdater<2>;
  template class GridOdePointIntegratorUpdater<3>;
}
