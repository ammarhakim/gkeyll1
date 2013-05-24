/**
 * @file	LcNodalHyperDiffusionUpdater.cpp
 *
 * @brief	Updater to evaluate (hyper)diffusion operators using nodal DG
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcNodalHyperDiffusionUpdater.h>

namespace Lucee
{
  template <> const char *NodalHyperDiffusionUpdater<1>::id = "HyperDiffusion1D";
  template <> const char *NodalHyperDiffusionUpdater<2>::id = "HyperDiffusion2D";
  template <> const char *NodalHyperDiffusionUpdater<3>::id = "HyperDiffusion3D";

  template <unsigned NDIM>
  NodalHyperDiffusionUpdater<NDIM>::NodalHyperDiffusionUpdater()
  {
  }

  template <unsigned NDIM>
  NodalHyperDiffusionUpdater<NDIM>::~NodalHyperDiffusionUpdater()
  {
  }

  template <unsigned NDIM>
  void
  NodalHyperDiffusionUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
  }

  template <unsigned NDIM>
  void
  NodalHyperDiffusionUpdater<NDIM>::initialize()
  {
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  NodalHyperDiffusionUpdater<NDIM>::update(double t)
  {
    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  NodalHyperDiffusionUpdater<NDIM>::declareTypes()
  {
  }

// instantiations
  template class NodalHyperDiffusionUpdater<1>;
  template class NodalHyperDiffusionUpdater<2>;
  template class NodalHyperDiffusionUpdater<3>;
}
