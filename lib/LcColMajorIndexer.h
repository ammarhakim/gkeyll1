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

namespace Lucee
{
/**
 * Base class provides key coefficient contruction functionality for
 * use in the indexing functions provided by the derived classes.
 */
  template <unsigned NDIM>
  struct ColMajorIndexerBase
  {
/**
 * Create a new indexer for mapping an N-dimensional index into a
 * linear index.
 *
 * @param start Starting index.
 * @param shape Shape of space.
 */
      ColMajorIndexerBase(unsigned shape[NDIM], int start[NDIM])
      {
        for (unsigned i=0; i<NDIM; ++i)
        {
          this->start[i] = start[i];
          this->shape[i] = shape[i];
        }

        ai[1] = 1;
        for (unsigned i=2; i<NDIM+1; ++i)
          ai[i] = ai[i-1]*shape[i-2];
        int sum = 0.0;
        for (unsigned i=1; i<NDIM+1; ++i)
          sum += ai[i]*start[i-1];
        ai[0] = -sum;
      }

/**
 * Return start index into space.
 *
 * @param i direction.
 * @return start index.
 */
      int getLower(unsigned i) const { return start[i]; }

/**
 * Return last index into space.
 *
 * @param i direction.
 * @return end index.
 */
      int getUpper(unsigned i) const { return start[i]+shape[i]; }

/**
 * Return linear index given N-dimensional index.
 *
 * @return Linear index.
 */
      int getGenIndex(int idx[NDIM]) const
      {
        int sum = ai[0]+idx[0];
        for (unsigned i=2; i<NDIM+1; ++i)
          sum += ai[i]*idx[i-1];
        return sum;
      }

    protected:
/** Coefficients for linear map */
      int ai[NDIM+1];

    private:
/** Start indices */
      int start[NDIM];
/** Shape of linear-space */
      unsigned shape[NDIM];
  };

/** 
 * Generic indexer class: empty except for base class methods.
 */
  template <unsigned NDIM>
  class ColMajorIndexer : public ColMajorIndexerBase<NDIM>
  {
    public:
/**
 * Create a new indexer.
 *
 * @param start Starting indices.
 * @param shape Shape of space.
 */
      ColMajorIndexer(unsigned shape[NDIM], int start[NDIM])
        : ColMajorIndexerBase<NDIM>(shape, start)
      {
      }
  };

/** One dimensional indexer */
  template <>
  class ColMajorIndexer<1> : public ColMajorIndexerBase<1>
  {
    public:
/**
 * Create a new indexer.
 *
 * @param start Starting indices.
 * @param shape Shape of space.
 */
      ColMajorIndexer(unsigned shape[1], int start[1])
        : ColMajorIndexerBase<1>(shape, start)
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
  class ColMajorIndexer<2> : public ColMajorIndexerBase<2>
  {
    public:
/**
 * Create a new indexer.
 *
 * @param start Starting indices.
 * @param shape Shape of space.
 */
      ColMajorIndexer(unsigned shape[2], int start[2])
        : ColMajorIndexerBase<2>(shape, start)
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
  class ColMajorIndexer<3> : public ColMajorIndexerBase<3>
  {
    public:
/**
 * Create a new indexer.
 *
 * @param start Starting indices.
 * @param shape Shape of space.
 */
      ColMajorIndexer(unsigned shape[3], int start[3])
        : ColMajorIndexerBase<3>(shape, start)
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
  class ColMajorIndexer<4> : public ColMajorIndexerBase<4>
  {
    public:
/**
 * Create a new indexer.
 *
 * @param start Starting indices.
 * @param shape Shape of space.
 */
      ColMajorIndexer(unsigned shape[4], int start[4])
        : ColMajorIndexerBase<4>(shape, start)
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

#endif // LC_COL_MAJOR_INDEXER_H
