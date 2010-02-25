/**
 * @file	LcColMajorSequencer.h
 *
 * @brief	Sequence over region using row-major order.
 *
 * @version	$Id: LcColMajorSequencer.h 278 2010-02-22 17:15:09Z a.hakim777 $
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

#ifndef LC_COL_MAJOR_SEQUENCER_H
#define LC_COL_MAJOR_SEQUENCER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcRegion.h>

namespace Lucee
{
  template <unsigned NDIM>
  class ColMajorSequencer
  {
    public:
/**
 * Create a new sequencer object to walk over a given region.
 *
 * @param rgn Region to sequence.
 */
      ColMajorSequencer(const Lucee::Region<NDIM, int>& rgn);

/**
 * Move sequencer by one step. Returns false there are no more steps
 * to take.
 *
 * @return true, if sequencing is not complete, false otherwise.
 */
      bool step();

/**
 * Reset the sequencer to point to the first location in region.
 */
      void reset();

/**
 * Return pointer to current indices in region.
 *
 * @return pointer to indices in region.
 */
      const int* getIndex() const { return idx; }

/**
 * Fill supplied array with current index in region.
 *
 * @param idx On output the current index in region.
 */
      void fillWithIndex(int idx[NDIM]) const;

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

  template <unsigned NDIM>
  ColMajorSequencer<NDIM>::ColMajorSequencer(const Lucee::Region<NDIM, int>& rgn)
    : rgn(rgn), isFirst(true)
  {
    isEmpty = rgn.getVolume() <= 0 ? true : false;
    for (unsigned i=0; i<NDIM; ++i)
      idx[i] = rgn.getLower(i);
  }

  template <unsigned NDIM>
  void
  ColMajorSequencer<NDIM>::fillWithIndex(int oidx[NDIM]) const
  {
    for (unsigned i=0; i<NDIM; ++i)
      oidx[i] = idx[i];
  }

  template <unsigned NDIM>
  bool
  ColMajorSequencer<NDIM>::step()
  {
    if (isEmpty)
    { // if region has zero volume, do nothing
      return false;
    }

    if (isFirst)
    { // first time around: indices already set in ctor
      isFirst = false;
      return true;
    }

    for (unsigned i=0; i<NDIM; ++i)
    {
      idx[i] += 1;
      if (idx[i] > rgn.getUpper(i)-1)
        idx[i] = rgn.getLower(i);
      else
        return true;
    }
    return false;
  }

  template <unsigned NDIM>
  void
  ColMajorSequencer<NDIM>::reset()
  {
    isFirst = true;
    for (unsigned i=0; i<NDIM; ++i)
      idx[i] = rgn.getLower(i);
  }
}

#endif // LC_COL_MAJOR_SEQUENCER_H
