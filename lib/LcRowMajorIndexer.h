/**
 * @file	LcRowMajorIndexer.h
 *
 * @brief	Class mapping an N-dimensional index into a linear index.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

#ifndef LC_ROW_MAJOR_INDEXER_H
#define LC_ROW_MAJOR_INDEXER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLinIndexer.h>
#include <LcRegion.h>
#include <LcRowMajorSequencer.h>

namespace Lucee
{
/**
 * Create a new indexer mapping coefficients.
 *
 * @param shape Shape of space.
 * @param start Start index.
 */
  template <unsigned NDIM>
  Lucee::FixedVector<NDIM+1, int>
  createRowMajorIndexer(const unsigned shape[NDIM], const int start[NDIM])
  {
    Lucee::FixedVector<NDIM+1, int> ai(0);
    ai[NDIM] = 1;
    for (unsigned i=NDIM-1; i>=1; --i)
      ai[i] = ai[i+1]*shape[i];

    int sum = 0;
    for (unsigned i=1; i<NDIM+1; ++i)
      sum += ai[i]*start[i-1];
    ai[0] = -sum;

    return ai;
  }

/**
 * Create a new indexer mapping coefficients.
 *
 * @param start Start index.
 * @param end Ending index.
 */
  template <unsigned NDIM>
  Lucee::FixedVector<NDIM+1, int>
  createRowMajorIndexer(const Lucee::Region<NDIM, int>& rgn)
  {
    Lucee::FixedVector<NDIM+1, int> ai(0);
    ai[NDIM] = 1;
    for (unsigned i=NDIM-1; i>=1; --i)
      ai[i] = ai[i+1]*rgn.getShape(i);

    int sum = 0;
    for (unsigned i=1; i<NDIM+1; ++i)
      sum += ai[i]*rgn.getLower(i-1);
    ai[0] = -sum;

    return ai;
  }

/**
 * Row-major indexer.
 */
  template <unsigned NDIM>
  class RowMajorIndexer : public Lucee::LinIndexer<NDIM>
  {
    public:
/** Type of sequencer that goes with this indexer */
      typedef RowMajorSequencer<NDIM> Sequencer;

/**
 * Create a new indexer for mapping an N-dimensional index into a
 * linear index.
 *
 * @param shape Shape of space.
 * @param start Starting index.
 */
      RowMajorIndexer(const unsigned shape[NDIM], const int start[NDIM])
        : Lucee::LinIndexer<NDIM>(shape, start, createRowMajorIndexer<NDIM>(shape, start))
      {
      }

/**
 * Create a new indexer over given N-dimensional region.
 *
 * @param rgn Region to index.
 */
      RowMajorIndexer(const Lucee::Region<NDIM, int>& rgn)
        : Lucee::LinIndexer<NDIM>(rgn, createRowMajorIndexer<NDIM>(rgn))
      {
      }

/**
 * Create a new indexer copying from input indexer.
 *
 * @param indexer Indexer to copy from.
 */
      RowMajorIndexer(const RowMajorIndexer<NDIM>& indexer)
        : Lucee::LinIndexer<NDIM>(indexer)
      {
      }

/**
 * Copy the values from the input indexer.
 *
 * @param indexer Indexer to copy from.
 * @return reference to this indexer.
 */
      RowMajorIndexer<NDIM>&
      operator=(const RowMajorIndexer<NDIM>& indexer)
      {
        if (&indexer == this)
          return *this;
        Lucee::LinIndexer<NDIM>::operator=(indexer);
        return *this;
      }
  };
}

#endif // LC_ROW_MAJOR_INDEXER_H
