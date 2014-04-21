/**
 * @file LcNonUniRectCartGrid.cpp
 *
 * @brief A grid allowing for non-uniform spacing, but otherwise rectangular
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcNonUniRectCartGrid.h>
#include <LcMathLib.h>
#include <LcStructGridField.h>

namespace Lucee
{
// set ids for grid creators
  template <> const char *NonUniRectCartGrid<1>::id = "NonUniformRectCart1D";
  template <> const char *NonUniRectCartGrid<2>::id = "NonUniformRectCart2D";
  template <> const char *NonUniRectCartGrid<3>::id = "NonUniformRectCart3D";

  template <unsigned NDIM>
  NonUniRectCartGrid<NDIM>::NonUniRectCartGrid()
  {
  }

  template <unsigned NDIM>
  NonUniRectCartGrid<NDIM>::~NonUniRectCartGrid()
  {
    for (unsigned d=0; d<NDIM; ++d)
    {
      delete vcoords[d];
      delete cellSize[d];
    }
  }

  template <unsigned NDIM>
  void
  NonUniRectCartGrid<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    StructuredGridBase<NDIM>::readInput(tbl);

// local region indexed by grid
    typename Lucee::Region<NDIM, int> localRgn = this->getLocalRegion();

// extend it to store approriate number of vertices
    Lucee::FixedVector<NDIM, int> lvExt(2), uvExt(3);
    typename Lucee::Region<NDIM, int> extVertRegion = localRgn.extend(&lvExt[0], &uvExt[0]);

// extend it to store approriate number of cells
    Lucee::FixedVector<NDIM, int> lcExt(2), ucExt(2);
    typename Lucee::Region<NDIM, int> extCellRegion = localRgn.extend(&lcExt[0], &ucExt[0]);

// allocate memory for storing coordinates of vertex
    for (unsigned d=0; d<NDIM; ++d)
    {
      int lower[1], upper[1];
      lower[0] = extVertRegion.getLower(d);
      upper[0] = extVertRegion.getUpper(d);
      Lucee::Region<1, int> rgn(lower, upper);
      vcoords[d] = new Array<1, double>(rgn);
    }

// get list of mapping functions
    Lucee::LuaTable mapTbl = tbl.getTable("mappings");

// allocate memory for storing cell size
    for (unsigned d=0; d<NDIM; ++d)
    {
      int lower[1], upper[1];
      lower[0] = extCellRegion.getLower(d);
      upper[0] = extCellRegion.getUpper(d);
      Lucee::Region<1, int> rgn(lower, upper);
      cellSize[d] = new Array<1, double>(rgn);
    }
  }

  template <unsigned NDIM>
  void
  NonUniRectCartGrid<NDIM>::getCentroid(double xc[3]) const
  {
  }

  template <unsigned NDIM>
  void
  NonUniRectCartGrid<NDIM>::getVertex(double xc[3]) const
  {
  }

  template <unsigned NDIM>
  double
  NonUniRectCartGrid<NDIM>::getVolume() const
  {
    return 0.0;
  }

  template <unsigned NDIM>
  double
  NonUniRectCartGrid<NDIM>::getSurfArea(unsigned dir) const
  {
    return 0.0;
  }

  template <unsigned NDIM>
  void
  NonUniRectCartGrid<NDIM>::getSurfCoordSys(unsigned dir, double norm[3],
    double tan1[3], double tan2[3]) const
  {
  }

  template <unsigned NDIM>
  TxIoNodeType
  NonUniRectCartGrid<NDIM>::writeToFile(TxIoBase& io, TxIoNodeType& node,
    const std::string& nm)
  {
// place-holder
    return node;
  }

  template <unsigned NDIM>
  void
  NonUniRectCartGrid<NDIM>::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
// call base class to register its methods
    Lucee::StructuredGridBase<NDIM>::appendLuaCallableMethods(lfm);
  }

// instantiations
  template class NonUniRectCartGrid<1>;
  template class NonUniRectCartGrid<2>;
  template class NonUniRectCartGrid<3>;
}
