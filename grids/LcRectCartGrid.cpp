/**
 * @file	LcRectCartGrid.cpp
 *
 * @brief	A rectangular cartesian grid.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcRectCartGrid.h>

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
  RectCartGrid<NDIM>::RectCartGrid(const Lucee::Region<NDIM, int>& localRgn,
    const Lucee::Region<NDIM, int>& globalRgn,
    const Lucee::Region<NDIM, double>& physBox) 
    : Lucee::StructuredGridBase<NDIM>(localRgn, globalRgn, physBox)
  {
    calcGeometry();
  }

  template <unsigned NDIM>
  void
  RectCartGrid<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    StructuredGridBase<NDIM>::readInput(tbl);
// compute stuff we need
    calcGeometry();
  }

  template <unsigned NDIM>
  void
  RectCartGrid<NDIM>::getCentroid(double xc[3]) const
  {
    for (unsigned i=0; i<NDIM; ++i)
      xc[i] = this->compSpace.getLower(i) +
        (this->currIdx[i]-this->globalRgn.getLower(i) + 0.5)*dx[i];
    for (unsigned i=NDIM; i<3; ++i)
      xc[i] = 0.0;
  }

  template <unsigned NDIM>
  void
  RectCartGrid<NDIM>::getVertex(double xc[3]) const
  {
    for (unsigned i=0; i<NDIM; ++i)
      xc[i] = this->compSpace.getLower(i) +
        (this->currIdx[i]-this->globalRgn.getLower(i))*dx[i];
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
  TxIoNodeType
  RectCartGrid<NDIM>::writeToFile(TxIoBase& io, TxIoNodeType& node,
    const std::string& nm)
  {
// create a group to write out grid data
    TxIoNodeType grdGrp = io.createGroup(node, nm);
    std::vector<double> lower(NDIM), upper(NDIM);
    std::vector<unsigned> numPhysCells(NDIM), start(NDIM);
    for (unsigned i=0; i<NDIM; ++i)
    {
      lower[i] = this->compSpace.getLower(i);
      upper[i] = this->compSpace.getUpper(i);
      start[i] = this->globalRgn.getLower(i);
      numPhysCells[i] = this->globalRgn.getShape(i);
    }
    io.writeAttribute(grdGrp, "vsType", "mesh");
    io.writeAttribute(grdGrp, "vsKind", "uniform");
    io.template 
      writeAttribute<unsigned>(grdGrp, "vsStartCell", start);
    io.template 
      writeAttribute<unsigned>(grdGrp, "vsNumCells", numPhysCells);
    io.template
      writeAttribute<double>(grdGrp, "vsLowerBounds", lower);
    io.template
      writeAttribute<double>(grdGrp, "vsUpperBounds", upper);

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

  template <unsigned NDIM>
  void
  RectCartGrid<NDIM>::calcGeometry()
  {
    for (unsigned i=0; i<3; ++i)
      dx[i] = 1.0;
    for (unsigned i=0; i<NDIM; ++i)
      dx[i] = this->compSpace.getShape(i)/this->globalRgn.getShape(i);
    cellVolume = 1.0;
    for (unsigned i=0; i<NDIM; ++i)
      cellVolume = cellVolume*dx[i];
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
