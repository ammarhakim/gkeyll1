/**
 * @file	LcGmvGridCreator.h
 *
 * @brief	Creator for grids defined in GMV format.
 *
 * @version	$Id$
 */

#ifndef LC_GMV_GRID_CREATOR_H
#define LC_GMV_GRID_CREATOR_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// std includes
#include <iostream>

// lucee includes
#include <LcUnstructGridCreator.h>

namespace Lucee
{
  template <typename REAL>
  class GmvGridCreator : public UnstructGridCreator<REAL>
  {
    public:
/**
 * Grid creator from a grid defined in GMV format.
 *
 * @param ndim Dimension of grid.
 * @param gmvStrm An input stream with GMV data.
 */
      GmvGridCreator(unsigned ndim, std::istream& gmvStrm);
  };
}

#endif // LC_GMV_GRID_CREATOR_H
