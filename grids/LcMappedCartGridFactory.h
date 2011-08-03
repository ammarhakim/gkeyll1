/**
 * @file	LcMappedCartGridFactory.h
 *
 * @brief	A factory to make mapped cartesian grids.
 */

#ifndef LC_MAPPED_CART_GRID_FACTORY_H
#define LC_MAPPED_CART_GRID_FACTORY_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLuaTable.h>

// std includes
#include <vector>

namespace Lucee
{
// forward declare MappedCartGrid
  template <unsigned NDIM> class MappedCartGrid;

/**
 * Factory to create rectangular cartesian grids.
 */
  template <unsigned NDIM>
  class MappedCartGridFactory
  {
    public:
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
      Lucee::MappedCartGrid<NDIM>* create();

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
