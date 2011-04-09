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
/**
 * Reader for meshes stored in GMV file. This is a very simple and
 * limited reader and does not support the full GMV format. It is only
 * meant for testing the grids classes and should *not* be used for
 * production use.
 */
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

    private:
/** Dimension of grid */
      unsigned ndim;

/**
 * Read grid from GMV format file.
 *
 * @param gmvStrm An input stream with GMV data.
 */
      void readFromGmvFile(std::istream& gmvStrm);
  };
}

#endif // LC_GMV_GRID_CREATOR_H
