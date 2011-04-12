/**
 * @file	LcCellElem.h
 *
 * @brief       Cell element class.
 *
 * @version	$Id$
 */

#ifndef LC_CELL_ELEM_H
#define LC_CELL_ELEM_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// std includes
#include <vector>

namespace Lucee
{
/**
 * A class representing a cell.
 */
  template <typename REAL>
  class CellElem
  {
    public:
/**
 * Fill with cell centroid coordinates.
 *
 * @param xv On output coordinates of cell centroid.
 */
      void fillWithCoordinates(REAL xv[3]) const;

/**
 * Get volume of cell.
 *
 * @return volume of cell.
 */
      REAL getMeasure() const;

    protected:
/**
 * Create cell given list of cell center coordinates and cell volumes.
 *
 * @param cc Cell center coordinates.
 * @param cv Cell volumes.
 */
      CellElem(const std::vector<REAL>& cc, const std::vector<REAL>& cv);

    private:
/** Reference to cell centroid coordinates */
      const std::vector<REAL>& cc;
/** Reference to cell volume */
      const std::vector<REAL>& cv;
/** Index of current cell */
      unsigned currCell;
  };
}

#endif // LC_CELL_ELEM_H
