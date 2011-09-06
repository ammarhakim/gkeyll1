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
// get input/output fields
    const Lucee::Field<NDIM, double>& inFld = this->getInp<Lucee::Field<NDIM, double> >(0);
    Lucee::Field<NDIM, double>& cdFld = this->getOut<Lucee::Field<NDIM, double> >(0);

    if (NDIM == 1)
      computeCentralDifference1D(inFld, cdFld);
    else if (NDIM == 2)
      computeCentralDifference2D(inFld, cdFld);
    else if (NDIM == 3)
      computeCentralDifference3D(inFld, cdFld);

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  RectSecondOrderCentralDiffUpdater<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void
  RectSecondOrderCentralDiffUpdater<NDIM>::computeCentralDifference1D(
    const Lucee::Field<NDIM, double>& inFld, Lucee::Field<NDIM, double>& cdFld)
  {
// get hold of grid
    const Lucee::RectCartGrid<1>& grid 
      = this->getGrid<Lucee::RectCartGrid<1> >();

    double dx2 = grid.getDx(0)*grid.getDx(0);

    for (int i=inFld.getLower(0); i<inFld.getUpper(0); ++i)
      for (unsigned k=0; k<inFld.getNumComponents(); ++k)
        cdFld(i,k) = (inFld(i-1,k) - 2.0*inFld(i,k) + inFld(i+1,k))/dx2;
  }

  template <unsigned NDIM>
  void
  RectSecondOrderCentralDiffUpdater<NDIM>::computeCentralDifference2D(
    const Lucee::Field<NDIM, double>& inFld, Lucee::Field<NDIM, double>& cdFld)
  {
// get hold of grid
    const Lucee::RectCartGrid<2>& grid 
      = this->getGrid<Lucee::RectCartGrid<2> >();

    double dx2 = grid.getDx(0)*grid.getDx(0);
    double dy2 = grid.getDx(1)*grid.getDx(1);

    for (int i=inFld.getLower(0); i<inFld.getUpper(0); ++i)
      for (int j=inFld.getLower(1); j<inFld.getUpper(1); ++j)
        for (unsigned k=0; k<inFld.getNumComponents(); ++k)
          cdFld(i,j,k) = (inFld(i-1,j,k) - 2.0*inFld(i,j,k) + inFld(i+1,j,k))/dx2
            + (inFld(i,j-1,k) - 2.0*inFld(i,j,k) + inFld(i,j+1,k))/dy2;
  }

  template <unsigned NDIM>
  void
  RectSecondOrderCentralDiffUpdater<NDIM>::computeCentralDifference3D(
    const Lucee::Field<NDIM, double>& inFld, Lucee::Field<NDIM, double>& cdFld)
  {
// get hold of grid
    const Lucee::RectCartGrid<2>& grid 
      = this->getGrid<Lucee::RectCartGrid<2> >();

    double dx2 = grid.getDx(0)*grid.getDx(0);
    double dy2 = grid.getDx(1)*grid.getDx(1);
    double dz2 = grid.getDx(2)*grid.getDx(2);

    for (int i=inFld.getLower(0); i<inFld.getUpper(0); ++i)
      for (int j=inFld.getLower(1); j<inFld.getUpper(1); ++j)
        for (int l=inFld.getLower(2); l<inFld.getUpper(1); ++l)
          for (unsigned k=0; k<inFld.getNumComponents(); ++k)
            cdFld(i,j,l,k) = (inFld(i-1,j,l,k) - 2.0*inFld(i,j,l,k) + inFld(i+1,j,l,k))/dx2
              + (inFld(i,j-1,l,k) - 2.0*inFld(i,j,l,k) + inFld(i,j+1,l,k))/dy2
              + (inFld(i,j,l-1,k) - 2.0*inFld(i,j,l,k) + inFld(i,j,l+1,k))/dz2;
  }

// instantiations
  template class RectSecondOrderCentralDiffUpdater<1>;
  template class RectSecondOrderCentralDiffUpdater<2>;
  template class RectSecondOrderCentralDiffUpdater<3>;
}
