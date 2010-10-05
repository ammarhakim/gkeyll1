/**
 * @file	LcRectCartGrid.cpp
 *
 * @brief	A rectangular cartesian grid.
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
#include <LcRectCartGrid.h>
#include <LcRectCartGridFactory.h>

namespace Lucee
{
// set ids for grid creators
  template <> const char *RectCartGrid<1>::id = "RectCart1D";
  template <> const char *RectCartGrid<2>::id = "RectCart2D";
  template <> const char *RectCartGrid<3>::id = "RectCart3D";

  template <unsigned NDIM>
  RectCartGrid<NDIM>::RectCartGrid()
  {
  }

  template <unsigned NDIM>
  RectCartGrid<NDIM>::RectCartGrid(const Lucee::Region<NDIM, int>& localBox,
    const Lucee::Region<NDIM, int>& globalBox,
    const Lucee::Region<NDIM, double>& physBox) 
    : Lucee::StructuredGridBase<NDIM>(localBox, globalBox, physBox)
  {
    for (unsigned i=0; i<3; ++i)
      dx[i] = 1.0;
    for (unsigned i=0; i<NDIM; ++i)
      dx[i] = physBox.getShape(i)/globalBox.getShape(i);
    cellVolume = 1.0;
    for (unsigned i=0; i<NDIM; ++i)
      cellVolume = cellVolume*dx[i];
  }

  template <unsigned NDIM>
  void
  RectCartGrid<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::RectCartGridFactory<NDIM> rgf;
    rgf.readInput(tbl);
    *this = *rgf.create();
  }

  template <unsigned NDIM>
  void
  RectCartGrid<NDIM>::getCentriod(double xc[3]) const
  {
    for (unsigned i=0; i<NDIM; ++i)
      xc[i] = this->compSpace.getLower(i) +
        (this->currIdx[i]-this->globalBox.getLower(i) + 0.5)*dx[i];
    for (unsigned i=NDIM; i<3; ++i)
      xc[i] = 0.0;
  }

  template <unsigned NDIM>
  void
  RectCartGrid<NDIM>::getVertex(double xc[3]) const
  {
    for (unsigned i=0; i<NDIM; ++i)
      xc[i] = (this->currIdx[i]-this->globalBox.getLower(i))*dx[i];
    for (unsigned i=NDIM; i<3; ++i)
      xc[i] = 0.0;
  }

  template <unsigned NDIM>
  double
  RectCartGrid<NDIM>::getVolume() const
  {
    return cellVolume;
  }

  template <unsigned NDIM>
  double
  RectCartGrid<NDIM>::getSurfArea(unsigned dir) const
  {
    if (dir==0)
      return dx[1]*dx[2];
    else if (dir==1)
      return dx[0]*dx[2];
    else if (dir==2)
      return dx[0]*dx[1];
    return 1.0;
  }

  template <unsigned NDIM>
  void
  RectCartGrid<NDIM>::getSurfCoordSys(unsigned dir, double norm[3],
    double tan1[3], double tan2[3]) const
  {
    if (dir==0)
    {
      norm[0] = 1.0;
      norm[1] = 0.0;
      norm[2] = 0.0;

      tan1[0] = 0.0;
      tan1[1] = 1.0;
      tan1[2] = 0.0;

      tan2[0] = 0.0;
      tan2[1] = 0.0;
      tan2[2] = 1.0;
    }
    else if (dir==1)
    {
      norm[0] = 0.0;
      norm[1] = 1.0;
      norm[2] = 0.0;

      tan1[0] = -1.0;
      tan1[1] = 0.0;
      tan1[2] = 0.0;

      tan2[0] = 0.0;
      tan2[1] = 0.0;
      tan2[2] = 1.0;
    }
    else if (dir==2)
    {
      norm[0] = 0.0;
      norm[1] = 0.0;
      norm[2] = 1.0;

      tan1[0] = 1.0;
      tan1[1] = 0.0;
      tan1[2] = 0.0;

      tan2[0] = 0.0;
      tan2[1] = 1.0;
      tan2[2] = 0.0;
    }
  }

  template <unsigned NDIM>
  Lucee::IoNodeType
  RectCartGrid<NDIM>::writeToFile(Lucee::IoBase& io, Lucee::IoNodeType& node,
    const std::string& nm)
  {
// create a group to write out grid data
    Lucee::IoNodeType grdGrp = io.createGroup(node, nm);
    std::vector<double> lower(NDIM), upper(NDIM);
    std::vector<unsigned> numPhysCells(NDIM), start(NDIM);
    for (unsigned i=0; i<NDIM; ++i)
    {
      lower[i] = this->compSpace.getLower(i);
      upper[i] = this->compSpace.getUpper(i);
      start[i] = this->globalBox.getLower(i);
      numPhysCells[i] = this->globalBox.getShape(i);
    }
    io.writeStrAttribute(grdGrp, "vsType", "mesh");
    io.writeStrAttribute(grdGrp, "vsKind", "uniform");
    io.template 
      writeVecAttribute<unsigned>(grdGrp, "vsStartCell", start);
    io.template 
      writeVecAttribute<unsigned>(grdGrp, "vsNumCells", numPhysCells);
    io.template
      writeVecAttribute<double>(grdGrp, "vsLowerBounds", lower);
    io.template
      writeVecAttribute<double>(grdGrp, "vsUpperBounds", upper);

    return grdGrp;
  }

  template <unsigned NDIM>
  void
  RectCartGrid<NDIM>::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
// call base class to register its methods
    Lucee::StructuredGridBase<NDIM>::appendLuaCallableMethods(lfm);
  }

  template <unsigned NDIM>
  RectCartGrid<NDIM>&
  RectCartGrid<NDIM>::operator=(const RectCartGrid<NDIM>& rg)
  {
    if (&rg == this)
      return *this;

    Lucee::StructuredGridBase<NDIM>::operator=(rg);
    for (unsigned i=0; i<3; ++i)
      dx[i] = rg.dx[i];
    cellVolume = rg.cellVolume;

    return *this;
  }

// instantiations
  template class RectCartGrid<1>;
  template class RectCartGrid<2>;
  template class RectCartGrid<3>;
  template class RectCartGrid<4>;
  template class RectCartGrid<5>;
  template class RectCartGrid<6>;
  template class RectCartGrid<7>;
}
