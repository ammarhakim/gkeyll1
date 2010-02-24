/**
 * @file	LcRowMajorSequencer.h
 *
 * @brief	Sequence over region using row-major order.
 *
 * @version	$Id: LcRowMajorSequencer.h 278 2010-02-22 17:15:09Z a.hakim777 $
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

#ifndef LC_ROW_MAJOR_SEQUENCER_H
#define LC_ROW_MAJOR_SEQUENCER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcRegion.h>

namespace Lucee
{
  template <unsigned NDIM>
  class RowMajorSequencer
  {
    public:
/**
 * Create a new sequencer object to walk over a given region.
 *
 * @param rgn Region to sequence.
 */
      RowMajorSequencer(const Lucee::Region<NDIM, int>& rgn);

/**
 * Move sequencer by one step. Returns false there are no more steps
 * to take.
 *
 * @return true, if sequencing is not complete, false otherwise.
 */
      bool step();

    private:
/** Region over which sequencer works */
      Lucee::Region<NDIM, int> rgn;
/** Flag to indicate if step() is called first time */
      bool isFirst;
/** Indices of current location */
      int idx[NDIM];
/** Is region empty */
      bool isEmpty;
  };
}

#endif // LC_ROW_MAJOR_SEQUENCER_H
