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
    : idxr(&Lucee::FixedVector<NDIM, unsigned>(1)[0], &Lucee::FixedVector<NDIM, int>(1)[0])
  {
  }

  template <unsigned NDIM>
  void
  NonUniRectCartGrid<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    StructuredGridBase<NDIM>::readInput(tbl);

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
    return ion;
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
