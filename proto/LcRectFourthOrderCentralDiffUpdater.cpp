/**
 * @file	LcRectFourthOrderCentralDiffUpdater.cpp
 *
 * @brief	Compute 2nd order central-differences on a rectangular grid.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcRectCartGrid.h>
#include <LcRectFourthOrderCentralDiffUpdater.h>
#include <LcStructGridField.h>

namespace Lucee
{
/** Class id: this is used by registration system */
  template <> const char *RectFourthOrderCentralDiffUpdater<1>::id = "RectFourthOrderCentralDiff1D";
  template <> const char *RectFourthOrderCentralDiffUpdater<2>::id = "RectFourthOrderCentralDiff2D";
  template <> const char *RectFourthOrderCentralDiffUpdater<3>::id = "RectFourthOrderCentralDiff3D";

  template <unsigned NDIM>
  RectFourthOrderCentralDiffUpdater<NDIM>::RectFourthOrderCentralDiffUpdater()
    : UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  RectFourthOrderCentralDiffUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  RectFourthOrderCentralDiffUpdater<NDIM>::update(double t)
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
  RectFourthOrderCentralDiffUpdater<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void
  RectFourthOrderCentralDiffUpdater<NDIM>::computeCentralDifference1D(
    const Lucee::Field<NDIM, double>& inFld, Lucee::Field<NDIM, double>& cdFld)
  {
// get hold of grid
    const Lucee::RectCartGrid<1>& grid 
      = this->getGrid<Lucee::RectCartGrid<1> >();

    double dx4 = grid.getDx(0)*grid.getDx(0)*grid.getDx(0)*grid.getDx(0);

    for (int i=inFld.getLower(0); i<inFld.getUpper(0); ++i)
      for (unsigned k=0; k<inFld.getNumComponents(); ++k)
        cdFld(i,k) = (inFld(i-2,k) - 4.0*inFld(i-1,k) + 6.0*inFld(i,k) - 4.0*inFld(i+1,k) + inFld(i+2,k))/dx4;
  }

  template <unsigned NDIM>
  void
  RectFourthOrderCentralDiffUpdater<NDIM>::computeCentralDifference2D(
    const Lucee::Field<NDIM, double>& inFld, Lucee::Field<NDIM, double>& cdFld)
  {
// get hold of grid
    const Lucee::RectCartGrid<2>& grid 
      = this->getGrid<Lucee::RectCartGrid<2> >();

    double dx4 = grid.getDx(0)*grid.getDx(0)*grid.getDx(0)*grid.getDx(0);
    double dy4 = grid.getDx(1)*grid.getDx(1)*grid.getDx(1)*grid.getDx(1);

    for (int i=inFld.getLower(0); i<inFld.getUpper(0); ++i)
      for (int j=inFld.getLower(1); j<inFld.getUpper(1); ++j)
        for (unsigned k=0; k<inFld.getNumComponents(); ++k)
          cdFld(i,j,k) =  (inFld(i-2,j,k) - 4.0*inFld(i-1,j,k) + 6.0*inFld(i,j,k) - 4.0*inFld(i+1,j,k) + inFld(i+2,j,k))/dx4 
            + (inFld(i,j-2,k) - 4.0*inFld(i,j-1,k) + 6.0*inFld(i,j,k) - 4.0*inFld(i,j+1,k) + inFld(i,j+2,k))/dy4;
  }

  template <unsigned NDIM>
  void
  RectFourthOrderCentralDiffUpdater<NDIM>::computeCentralDifference3D(
    const Lucee::Field<NDIM, double>& inFld, Lucee::Field<NDIM, double>& cdFld)
  {
// get hold of grid
    const Lucee::RectCartGrid<3>& grid 
      = this->getGrid<Lucee::RectCartGrid<3> >();

    double dx4 = grid.getDx(0)*grid.getDx(0)*grid.getDx(0)*grid.getDx(0);
    double dy4 = grid.getDx(1)*grid.getDx(1)*grid.getDx(1)*grid.getDx(1);
    double dz4 = grid.getDx(2)*grid.getDx(2)*grid.getDx(2)*grid.getDx(2);

    for (int i=inFld.getLower(0); i<inFld.getUpper(0); ++i)
      for (int j=inFld.getLower(1); j<inFld.getUpper(1); ++j)
        for (int l=inFld.getLower(2); l<inFld.getUpper(2); ++l)
          for (unsigned k=0; k<inFld.getNumComponents(); ++k)
            cdFld(i,j,l,k) =(inFld(i-2,j,l,k) - 4.0*inFld(i-1,j,l,k) + 6.0*inFld(i,j,l,k) - 4.0*inFld(i+1,j,l,k) + inFld(i+2,j,l,k))/dx4 
              + (inFld(i,j-2,l,k) - 4.0*inFld(i,j-1,l,k) + 6.0*inFld(i,j,l,k) - 4.0*inFld(i,j+1,l,k) + inFld(i,j+2,l,k))/dy4 
              + (inFld(i,j,l-2,k) - 4.0*inFld(i,j,l-1,k) + 6.0*inFld(i,j,l,k) - 4.0*inFld(i,j,l+1,k) + inFld(i,j,l+2,k))/dz4;
  }

// instantiations
  template class RectFourthOrderCentralDiffUpdater<1>;
  template class RectFourthOrderCentralDiffUpdater<2>;
  template class RectFourthOrderCentralDiffUpdater<3>;
}
