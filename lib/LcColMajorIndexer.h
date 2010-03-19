/**
 * @file	LcColMajorIndexer.h
 *
 * @brief	Class mapping an N-dimensional index into a linear index.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

#ifndef LC_COL_MAJOR_INDEXER_H
#define LC_COL_MAJOR_INDEXER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcColMajorSequencer.h>
#include <LcLinIndexer.h>
#include <LcRegion.h>

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
  createColMajorIndexer(const unsigned shape[NDIM], const int start[NDIM])
  {
    Lucee::FixedVector<NDIM+1, int> ai(0);
    ai[1] = 1;
    for (unsigned i=2; i<NDIM+1; ++i)
      ai[i] = ai[i-1]*shape[i-2];
    int sum = 0.0;
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
  createColMajorIndexer(const Lucee::Region<NDIM, int>& rgn)
  {
    Lucee::FixedVector<NDIM+1, int> ai(0);
    ai[1] = 1;
    for (unsigned i=2; i<NDIM+1; ++i)
      ai[i] = ai[i-1]*rgn.getShape(i-2);
    int sum = 0.0;
    for (unsigned i=1; i<NDIM+1; ++i)
      sum += ai[i]*rgn.getLower(i-1);
    ai[0] = -sum;

    return ai;
  }

/**
 * Col-major indexer.
 */
  template <unsigned NDIM>
  class ColMajorIndexer : public Lucee::LinIndexer<NDIM>
  {
    public:
/** Type of sequencer that goes with this indexer */
      typedef ColMajorSequencer<NDIM> Sequencer;

/**
 * Create a new indexer for mapping an N-dimensional index into a
 * linear index.
 *
 * @param shape Shape of space.
 * @param start Starting index.
 */
      ColMajorIndexer(const unsigned shape[NDIM], const int start[NDIM])
        : Lucee::LinIndexer<NDIM>(shape, start, createColMajorIndexer<NDIM>(shape, start))
      {
      }

/**
 * Create a new indexer over given N-dimensional region.
 *
 * @param rgn Region to index.
 */
      ColMajorIndexer(const Lucee::Region<NDIM, int>& rgn)
        : Lucee::LinIndexer<NDIM>(rgn, createColMajorIndexer<NDIM>(rgn))
      {
      }

/**
 * Create a new indexer copying from input indexer.
 *
 * @param indexer Indexer to copy from.
 */
      ColMajorIndexer(const ColMajorIndexer<NDIM>& indexer)
        : Lucee::LinIndexer<NDIM>(indexer)
      {
      }

/**
 * Create a new indexer copying from input indexer.
 *
 * @param indexer Indexer to copy from.
 */
      ColMajorIndexer(const LinIndexer<NDIM>& indexer)
        : Lucee::LinIndexer<NDIM>(indexer)
      {
      }

/**
 * Copy the values from the input indexer.
 *
 * @param indexer Indexer to copy from.
 * @return reference to this indexer.
 */
      ColMajorIndexer<NDIM>&
      operator=(const ColMajorIndexer<NDIM>& indexer)
      {
        if (&indexer == this)
          return *this;
        Lucee::LinIndexer<NDIM>::operator=(indexer);
        return *this;
      }

/**
 * Deflate the index into a lower-dimensional indexer. The returned
 * array shares the data with this indexer.
 *
 * @param defDims Dimensions to remove from indexer.
 * @param defDimsIdx Index along the removed dimensions.
 * @return deflated indexer.
 */
      template <unsigned RDIM>
      ColMajorIndexer<RDIM>
      deflate(const unsigned defDims[NDIM-RDIM], const int defDimsIdx[NDIM-RDIM]) const
      {
        return this->template deflateLin<RDIM>(defDims, defDimsIdx);
      }
  };
}

#endif // LC_COL_MAJOR_INDEXER_H
