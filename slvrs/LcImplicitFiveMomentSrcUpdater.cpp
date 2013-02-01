/**
 * @file	LcImplicitFiveMomentSrcUpdater.cpp
 *
 * @brief	Implicit updater for 5-moment source terms
 */

// lucee includes
#include <LcImplicitFiveMomentSrcUpdater.h>

namespace Lucee
{
// set ids for module system
  template <> const char *ImplicitFiveMomentSrcUpdater<1>::id = "ImplicitFiveMomentSrc1D";
  template <> const char *ImplicitFiveMomentSrcUpdater<2>::id = "ImplicitFiveMomentSrc2D";
  template <> const char *ImplicitFiveMomentSrcUpdater<3>::id = "ImplicitFiveMomentSrc3D";

  template <unsigned NDIM>
  void
  ImplicitFiveMomentSrcUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    UpdaterIfc::readInput(tbl);
  }

  template <unsigned NDIM>
  void
  ImplicitFiveMomentSrcUpdater<NDIM>::initialize()
  {
// call base class method
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  ImplicitFiveMomentSrcUpdater<NDIM>::update(double t)
  {
    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  ImplicitFiveMomentSrcUpdater<NDIM>::declareTypes()
  {
    this->setLastOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class ImplicitFiveMomentSrcUpdater<1>;
  template class ImplicitFiveMomentSrcUpdater<2>;
  template class ImplicitFiveMomentSrcUpdater<3>;
}

