/**
 * @file   LcGridOdeIntegrator.cpp
 *
 * @brief   Base class for ODE integrator over complete grid.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGridOdeIntegrator.h>

namespace Lucee
{
  template <unsigned NDIM>
  GridOdeIntegrator<NDIM>::GridOdeIntegrator(
    const Lucee::StructuredGridBase<NDIM>& grid, unsigned numInp)
    : gridPtr(&grid), numInp(numInp), inpFlds(numInp)
  {
  }

  template <unsigned NDIM>
  void
  GridOdeIntegrator<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
  }

  template <unsigned NDIM>
  const Lucee::Field<NDIM, double>&
  GridOdeIntegrator<NDIM>::getIn(unsigned loc) const
  {
    if (loc >= inpFlds.size())
    {
      Lucee::Except lce("GridOdeIntegrator::getIn: Incorrect index ");
      lce << loc << " specified. Index should be less than " << inpFlds.size();
      throw lce;
    }
    return *inpFlds[loc];
  }

  template <unsigned NDIM>
  void
  GridOdeIntegrator<NDIM>::setIn(unsigned loc, const Lucee::Field<NDIM, double>& in)
  {
    if (loc >= inpFlds.size())
    {
      Lucee::Except lce("GridOdeIntegrator::setIn: Incorrect index ");
      lce << loc << " specified. Index should be less than " << inpFlds.size();
      throw lce;
    }
    inpFlds[loc] = &in;
  }

// instantiations
  template class GridOdeIntegrator<1>;
  template class GridOdeIntegrator<2>;
  template class GridOdeIntegrator<3>;
}
