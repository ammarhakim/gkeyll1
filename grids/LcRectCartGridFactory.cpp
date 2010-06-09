/**
 * @file	LcRectCartGridFactory.cpp
 *
 * @brief	A facotry to make rectangular cartesian grids.
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
#include <LcExcept.h>
#include <LcRectCartGridFactory.h>
#include <LcRectCartGrid.h>

namespace Lucee
{
// set ids for grid creators
  template <> const char *RectCartGridFactory<1>::id = "RectCart1D";
  template <> const char *RectCartGridFactory<2>::id = "RectCart2D";
  template <> const char *RectCartGridFactory<3>::id = "RectCart3D";
  
  template <unsigned NDIM>
  void
  RectCartGridFactory<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    cells = tbl.getNumVec("cells");
    if (cells.size() != NDIM)
    {
      Lucee::Except lce("RectCartGridFactory::readInput: 'cells' should have exactly ");
      lce << NDIM << " elements. Instead has " << cells.size() << std::endl;
    }

    lower = tbl.getNumVec("lower");
    if (lower.size() != NDIM)
    {
      Lucee::Except lce("RectCartGridFactory::readInput: 'lower' should have exactly ");
      lce << NDIM << " elements. Instead has " << lower.size() << std::endl;
    }

    upper = tbl.getNumVec("upper");
    if (upper.size() != NDIM)
    {
      Lucee::Except lce("RectCartGridFactory::readInput: 'upper' should have exactly ");
      lce << NDIM << " elements. Instead has " << upper.size() << std::endl;
    }
  }

  template <unsigned NDIM>
  GridIfc*
  RectCartGridFactory<NDIM>::create(const Lucee::SolverIfc& solver)
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
    return new RectCartGrid<NDIM>(lrgn, grgn, dom); 
  }

// instantiations
  template class RectCartGridFactory<1>;
  template class RectCartGridFactory<2>;
  template class RectCartGridFactory<3>;
}
