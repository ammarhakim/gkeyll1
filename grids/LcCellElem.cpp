/**
 * @file	LcCellElem.cpp
 *
 * @brief       Cell element class.
 *
 * @version	$Id$
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCellElem.h>

namespace Lucee
{
  template <typename REAL>
  void
  CellElem<REAL>::fillWithCoordinates(REAL xv[3]) const
  {
    for (unsigned i=0; i<3; ++i)
      xv[i] = cc[3*currCell+i];
  }

  template <typename REAL>
  REAL
  CellElem<REAL>::getMeasure() const
  {
    return cv[currCell];
  }

  template <typename REAL>
  CellElem<REAL>::CellElem(const std::vector<REAL>& cc, const std::vector<REAL>& cv)
    : cc(cc), cv(cv), currCell(0)
  {
  }

// instantiations
  template class CellElem<float>;
  template class CellElem<double>;
}
