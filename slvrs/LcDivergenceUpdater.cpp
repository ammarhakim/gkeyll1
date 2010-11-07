/**
 * @file	LcDivergenceUpdater.cpp
 *
 * @brief	Compute divergence on rectangular grids.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDivergenceUpdater.h>
#include <LcField.h>

namespace Lucee
{
  template <> const char *DivergenceUpdater<1>::id = "Divergence1D";
  template <> const char *DivergenceUpdater<2>::id = "Divergence2D";
  template <> const char *DivergenceUpdater<3>::id = "Divergence3D";

  template <unsigned NDIM>
  DivergenceUpdater<NDIM>::DivergenceUpdater()
    : UpdaterIfc() 
  {
  }

  template <unsigned NDIM>
  void
  DivergenceUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  DivergenceUpdater<NDIM>::update(double t)
  {
    

    return UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  DivergenceUpdater<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class DivergenceUpdater<1>;
  template class DivergenceUpdater<2>;
  template class DivergenceUpdater<3>;
}
