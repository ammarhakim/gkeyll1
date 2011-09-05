/**
 * @file	LcRectSecondOrderCentralDiffUpdater.cpp
 *
 * @brief	Compute 2nd order central-differences on a rectangular grid.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcRectCartGrid.h>
#include <LcRectSecondOrderCentralDiffUpdater.h>
#include <LcStructGridField.h>

namespace Lucee
{
/** Class id: this is used by registration system */
  template <> const char *RectSecondOrderCentralDiffUpdater<1>::id = "RectSecondOrderCentralDiff1D";
  template <> const char *RectSecondOrderCentralDiffUpdater<2>::id = "RectSecondOrderCentralDiff2D";
  template <> const char *RectSecondOrderCentralDiffUpdater<3>::id = "RectSecondOrderCentralDiff3D";

  template <unsigned NDIM>
  RectSecondOrderCentralDiffUpdater<NDIM>::RectSecondOrderCentralDiffUpdater()
    : UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  RectSecondOrderCentralDiffUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  RectSecondOrderCentralDiffUpdater<NDIM>::update(double t)
  {
// get hold of grid
    const Lucee::RectCartGrid<NDIM>& grid 
      = this->getGrid<Lucee::RectCartGrid<NDIM> >();
// get input/output fields
    const Lucee::Field<NDIM, double>& infld = this->getInp<Lucee::Field<NDIM, double> >(0);
    Lucee::Field<NDIM, double>& cdFld = this->getOut<Lucee::Field<NDIM, double> >(0);


    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  RectSecondOrderCentralDiffUpdater<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class RectSecondOrderCentralDiffUpdater<1>;
  template class RectSecondOrderCentralDiffUpdater<2>;
  template class RectSecondOrderCentralDiffUpdater<3>;
}
