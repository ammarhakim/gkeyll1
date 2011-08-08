/**
 * @file	LcMappedCartGridFactory.cpp
 *
 * @brief	A facotry to make rectangular cartesian grids.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcExcept.h>
#include <LcMappedCartGridFactory.h>
#include <LcMappedCartGrid.h>

namespace Lucee
{
  template <unsigned NDIM>
  void
  MappedCartGridFactory<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    cells = tbl.getNumVec("cells");
    if (cells.size() != NDIM)
    {
      Lucee::Except lce("MappedCartGridFactory::readInput: 'cells' should have exactly ");
      lce << NDIM << " elements. Instead has " << cells.size() << std::endl;
      throw lce;
    }

    lower = tbl.getNumVec("lower");
    if (lower.size() != NDIM)
    {
      Lucee::Except lce("MappedCartGridFactory::readInput: 'lower' should have exactly ");
      lce << NDIM << " elements. Instead has " << lower.size() << std::endl;
      throw lce;
    }

    upper = tbl.getNumVec("upper");
    if (upper.size() != NDIM)
    {
      Lucee::Except lce("MappedCartGridFactory::readInput: 'upper' should have exactly ");
      lce << NDIM << " elements. Instead has " << upper.size() << std::endl;
      throw lce;
    }
  }

  template <unsigned NDIM>
  Lucee::MappedCartGrid<NDIM>*
  MappedCartGridFactory<NDIM>::create()
  {
    int ilo[NDIM], iup[NDIM];
    double xlo[NDIM], xup[NDIM];
    for (unsigned i=0; i<NDIM; ++i)
    {
      ilo[i] = 0;
      iup[i] = cells[i];
      xlo[i] = lower[i];
      xup[i] = upper[i];
    }
    Lucee::Region<NDIM, int> lrgn(ilo, iup);
    Lucee::Region<NDIM, int> grgn(ilo, iup);
    Lucee::Region<NDIM, double> dom(xlo, xup);
// ghost cells required to store geometry in extended region
    int lg[NDIM], ug[NDIM];
    for (unsigned i=0; i<NDIM; ++i)
    {
      lg[i] = 1;
      ug[i] = 2;
    }
    return new Lucee::MappedCartGrid<NDIM>(lrgn, lrgn.extend(lg, ug),
      grgn, dom);
// THE NODAL COORDINATES NEED TO BE COMPUTED HERE AND SENT TO CONSTRUCTOR
  }

// instantiations
  template class MappedCartGridFactory<1>;
  template class MappedCartGridFactory<2>;
  template class MappedCartGridFactory<3>;
}
