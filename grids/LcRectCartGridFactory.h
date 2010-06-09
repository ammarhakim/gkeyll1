/**
 * @file	LcRectCartGridFactory.h
 *
 * @brief	A facotry to make rectangular cartesian grids.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_RECT_CART_GRID_FACTORY_H
#define LC_RECT_CART_GRID_FACTORY_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGenericFactory.h>
#include <LcGridBase.h>
#include <LcRegion.h>

// std includes
#include <vector>

namespace Lucee
{
  template <unsigned NDIM>
  class RectCartGridFactory : public Lucee::GenericFactory<Lucee::GridBase>
  {
    public:
/** Class id: this is used by the registration system */
      static const char *id;

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      void readInput(Lucee::LuaTable& tbl);

/**
 * Create a new rectangular cartesian grid.
 *
 * @return pointer to new grid.
 */
      Lucee::GridBase* create(const Lucee::SolverIfc& solver);

    private:
/** Number of cells in domain */
      std::vector<double> cells;
/** Lower coordinates of space */
      std::vector<double> lower;
/** Upper coordinates of space */
      std::vector<double> upper;
  };
}

#endif // LC_RECT_CART_GRID_FACTORY_H
