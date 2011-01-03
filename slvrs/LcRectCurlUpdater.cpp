/**
 * @file	LcRectCurlUpdater.cpp
 *
 * @brief	Compute curl on rectangular grids.
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
#include <LcField.h>
#include <LcRectCurlUpdater.h>

namespace Lucee
{
  template <> const char *RectCurlUpdater<1>::id = "Curl1D";
  template <> const char *RectCurlUpdater<2>::id = "Curl2D";
  template <> const char *RectCurlUpdater<3>::id = "Curl3D";

  template <unsigned NDIM>
  RectCurlUpdater<NDIM>::RectCurlUpdater()
    : UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  RectCurlUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// read multiplication factor
    alpha = tbl.getNumber("alpha");
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  RectCurlUpdater<NDIM>::update(double t)
  {

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  RectCurlUpdater<NDIM>::declareTypes()
  {
// two input fields
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
// one output field
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class RectCurlUpdater<1>;
  template class RectCurlUpdater<2>;
  template class RectCurlUpdater<3>;
}
