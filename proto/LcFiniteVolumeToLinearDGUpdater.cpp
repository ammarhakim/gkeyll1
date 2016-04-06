/**
 * @file	LcFiniteVolumeToLinearDGUpdater.cpp
 *
 * @brief	Updater to integrate nodal DG/CG field over complete domain
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGlobals.h>
#include <LcFiniteVolumeToLinearDGUpdater.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  template <> const char *FiniteVolumeToLinearDGUpdater<1>::id = "FiniteVolumeToLinearDG1D";
  template <> const char *FiniteVolumeToLinearDGUpdater<2>::id = "FiniteVolumeToLinearDG2D";
  template <> const char *FiniteVolumeToLinearDGUpdater<3>::id = "FiniteVolumeToLinearDG3D";

  template <unsigned NDIM>
  FiniteVolumeToLinearDGUpdater<NDIM>::FiniteVolumeToLinearDGUpdater()
    : Lucee::UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  FiniteVolumeToLinearDGUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  FiniteVolumeToLinearDGUpdater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    const Lucee::Field<NDIM, double>& fvFld = this->getInp<Lucee::Field<NDIM, double> >(0);
    Lucee::Field<NDIM, double>& dgFld = this->getOut<Lucee::Field<NDIM, double> >(0);

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  FiniteVolumeToLinearDGUpdater<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class Lucee::FiniteVolumeToLinearDGUpdater<1>;
  template class Lucee::FiniteVolumeToLinearDGUpdater<2>;
  template class Lucee::FiniteVolumeToLinearDGUpdater<3>;
}
