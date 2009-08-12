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

namespace Lucee
{
/**
 * Base class provides key coefficient contruction functionality for
 * use in the indexing functions provided by the derived classes.
 */
  template <unsigned NDIM>
  struct RowMajorIndexerBase
  {
/**
 * Create a new indexer for mapping an N-dimensional index into a
 * linear index.
 *
 * @param start Starting index.
 * @param shape Shape of space.
 */
      RowMajorIndexerBase(int start[NDIM], unsigned shape[NDIM])
      {
        a[1] = 1;
        for (unsigned i=2; i<NDIM+1; ++i)
          ai[i] = ai[i-1]*shape[i-2];
        int sum = 0.0;
        for (unsigned i=1; i<NDIM+1; ++i)
          sum += a[i]*start[i-1];
        a[0] = -sum;
      }

/**
 * Return linear index given N-dimensional index.
 *
 * @return Linear index.
 */
      int getIndex(int idx[NDIM]) const 
      {
        int sum = ai[0]+idx[0];
        for (unsigned i=2; i<NDIM+1; ++i)
          sum += ai[i]*idx[i-1];
        return sum;
      }

/** Coefficients for linear map */
      int ai[NDIM+1];
  };

/** 
 * Generic indexer class: empty except for base class methods.
 */
  template <unsigned NDIM>
  class RowMajorIndexer : public RowMajorIndexerBase<NDIM>
  {
    public:
/**
 * Create a new indexer.
 *
 * @param start Starting indices.
 * @param shape Shape of space.
 */
      RowMajorIndexer(int start[NDIM], unsigned shape[NDIM])
        : RowMajorIndexerBase<NDIM>(start, shape)
      {
      }
  };

/** One dimensional indexer */
  template <>
  class RowMajorIndexer<1> : public RowMajorIndexerBase<1>
  {
    public:
/**
 * Create a new indexer.
 *
 * @param start Starting indices.
 * @param shape Shape of space.
 */
      RowMajorIndexer(int start[1], unsigned shape[1])
        : RowMajorIndexerBase<1>(start, shape)
      {
      }

/**
 * Map 1D index to a linear index.
 *
 * @param i Index location.
 * @return Index of (i) into linear space.
 */
      int getIndex(int i) const {
        return ai[0]+i;
      }
  };

/** Two dimensional indexer */
  template <>
  class RowMajorIndexer<2> : public RowMajorIndexerBase<2>
  {
    public:
/**
 * Create a new indexer.
 *
 * @param start Starting indices.
 * @param shape Shape of space.
 */
      RowMajorIndexer(int start[2], unsigned shape[2])
        : RowMajorIndexerBase<2>(start, shape)
      {
      }

/**
 * Map 2D index to a linear index.
 *
 * @param i Index location.
 * @param j Index location.
 * @return Index of (i,j) into linear space.
 */
      int getIndex(int i, int j) const {
        return ai[0]+i+ai[2]*j;
      }
  };

/** Three dimensional indexer */
  template <>
  class RowMajorIndexer<3> : public RowMajorIndexerBase<3>
  {
    public:
/**
 * Create a new indexer.
 *
 * @param start Starting indices.
 * @param shape Shape of space.
 */
      RowMajorIndexer(int start[3], unsigned shape[3])
        : RowMajorIndexerBase<3>(start, shape)
      {
      }

/**
 * Map 3D index to a linear index.
 *
 * @param i Index location.
 * @param j Index location.
 * @param k Index location.
 * @return Index of (i,j,k) into linear space.
 */
      int getIndex(int i, int j, int k) const {
        return ai[0]+i+ai[2]*j+ai[3]*k;
      }
  };

/** Four dimensional indexer */
  template <>
  class RowMajorIndexer<4> : public RowMajorIndexerBase<4>
  {
    public:
/**
 * Create a new indexer.
 *
 * @param start Starting indices.
 * @param shape Shape of space.
 */
      RowMajorIndexer(int start[4], unsigned shape[4])
        : RowMajorIndexerBase<4>(start, shape)
      {
      }

/**
 * Map 4D index to a linear index.
 *
 * @param i Index location.
 * @param j Index location.
 * @param k Index location.
 * @param l Index location.
 * @return Index of (i,j,k,l) into linear space.
 */
      int getIndex(int i, int j, int k, int l) const {
        return ai[0]+i+ai[2]*j+ai[3]*k+ai[4]*l;
      }
  };
}

#endif // LC_ROW_MAJOR_INDEXER_H
