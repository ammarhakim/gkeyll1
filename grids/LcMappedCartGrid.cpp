/**
 * @file LcMappedCartGrid.cpp
 *
 * @brief A logically rectangular grid, but non-rectangular in physical space.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcMappedCartGrid.h>

namespace Lucee
{
// set ids for grid creators
  template <> const char *MappedCartGrid<1>::id = "MappedCart1D";
  template <> const char *MappedCartGrid<2>::id = "MappedCart2D";
  template <> const char *MappedCartGrid<3>::id = "MappedCart3D";

  template <unsigned NDIM>
  MappedCartGrid<NDIM>::MappedCartGrid()
  {
  }

  template <unsigned NDIM>
  MappedCartGrid<NDIM>::MappedCartGrid(const Lucee::Region<NDIM, int>& localBox,
    const Lucee::Region<NDIM, int>& localExtBox,
    const Lucee::Region<NDIM, int>& globalBox,
    const Lucee::Region<NDIM, double>& compSpace) 
    : Lucee::StructuredGridBase<NDIM>(localBox, globalBox, compSpace),
      localExtBox(localExtBox)
  {
  }

  template <unsigned NDIM>
  void
  MappedCartGrid<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
  }

  template <unsigned NDIM>
  void
  MappedCartGrid<NDIM>::getCentriod(double xc[3]) const
  {
  }

  template <unsigned NDIM>
  void
  MappedCartGrid<NDIM>::getVertex(double xc[3]) const
  {
  }

  template <unsigned NDIM>
  double
  MappedCartGrid<NDIM>::getVolume() const
  {
    return 0.0;
  }

  template <unsigned NDIM>
  double
  MappedCartGrid<NDIM>::getSurfArea(unsigned dir) const
  {
    return 0.0;
  }

  template <unsigned NDIM>
  void
  MappedCartGrid<NDIM>::getSurfCoordSys(unsigned dir, double norm[3],
    double tan1[3], double tan2[3]) const
  {
  }

  template <unsigned NDIM>
  Lucee::IoNodeType
  MappedCartGrid<NDIM>::writeToFile(Lucee::IoBase& io, Lucee::IoNodeType& node,
    const std::string& nm)
  {
    return node;
  }

  template <unsigned NDIM>
  void
  MappedCartGrid<NDIM>::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
// call base class to register its methods
    Lucee::StructuredGridBase<NDIM>::appendLuaCallableMethods(lfm);
  }

// instantiations
  template class MappedCartGrid<1>;
  template class MappedCartGrid<2>;
  template class MappedCartGrid<3>;
}
