/**
 * @file	LcArrayIndexer.h
 *
 * @brief	Class mapping an N-dimensional index into a linear index.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

#ifndef LC_ARRAY_INDEXER_H
#define LC_ARRAY_INDEXER_H

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
  struct ArrayIndexerBase
  {
/**
 * Create a new indexer for mapping an N-dimensional index into a
 * linear index.
 *
 * @param start Starting index.
 * @param shape Shape of space.
 */
      ArrayIndexerBase(int start[NDIM], unsigned shape[NDIM])
      {
      }

/**
 * Return linear index given N-dimensional index.
 *
 * @return Linear index.
 */
      int getIndex(int idx[NDIM]) const 
      {
        int sum = ai[0];
        for (unsigned i=0; i<NDIM; ++i)
          sum += ai[i]*idx[i-1];
        return sum;
      }

/** Coefficients used for linear map */
      int ai[NDIM+1];
  };

/** 
 * Generic indexer class: empty except for base class methods.
 */
  template <unsigned NDIM>
  class ArrayIndexer : public ArrayIndexerBase<NDIM>
  {
    public:
/**
 * Create a new indexer.
 *
 * @param start Starting indices.
 * @param shape Shape of space.
 */
      ArrayIndexer(int start[NDIM], unsigned shape[NDIM])
        : ArrayIndexerBase<NDIM>(start, shape)
      {
      }
  };

/** One dimensional indexer */
  template <>
  class ArrayIndexer<1> : public ArrayIndexerBase<1>
  {
    public:
/**
 * Create a new indexer.
 *
 * @param start Starting indices.
 * @param shape Shape of space.
 */
      ArrayIndexer(int start[1], unsigned shape[1])
        : ArrayIndexerBase<1>(start, shape)
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
  class ArrayIndexer<2> : public ArrayIndexerBase<2>
  {
    public:
/**
 * Create a new indexer.
 *
 * @param start Starting indices.
 * @param shape Shape of space.
 */
      ArrayIndexer(int start[2], unsigned shape[2])
        : ArrayIndexerBase<2>(start, shape)
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
  class ArrayIndexer<3> : public ArrayIndexerBase<3>
  {
    public:
/**
 * Create a new indexer.
 *
 * @param start Starting indices.
 * @param shape Shape of space.
 */
      ArrayIndexer(int start[3], unsigned shape[3])
        : ArrayIndexerBase<3>(start, shape)
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
  class ArrayIndexer<4> : public ArrayIndexerBase<4>
  {
    public:
/**
 * Create a new indexer.
 *
 * @param start Starting indices.
 * @param shape Shape of space.
 */
      ArrayIndexer(int start[4], unsigned shape[4])
        : ArrayIndexerBase<4>(start, shape)
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

#endif // LC_ARRAY_INDEXER_H
