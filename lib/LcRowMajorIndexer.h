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
#include <LcRegion.h>

namespace Lucee
{
/**
 * Base class provides key coefficient contruction functionality for
 * use in the indexing functions provided by the derived classes.
 */
  template <unsigned NDIM>
  class RowMajorIndexerBase
  {
    public:
/**
 * Create a new indexer for mapping an N-dimensional index into a
 * linear index.
 *
 * @param shp Shape of space.
 * @param sta Starting index.
 */
      RowMajorIndexerBase(const unsigned shp[NDIM], const int sta[NDIM])
      {
        createRowMajorIndexer(shp, sta);
      }

/**
 * Create a new indexer over given N-dimensional region.
 *
 * @param rgn Region to index.
 */
      RowMajorIndexerBase(const Region<NDIM, int>& rgn)
      {
        unsigned shp[NDIM];
        int sta[NDIM];
        for (unsigned i=0; i<NDIM; ++i)
        {
          shp[i] = rgn.getShape(i);
          sta[i] = rgn.getLower(i);
        }
        createRowMajorIndexer(shp, sta);
      }

/**
 * Create a new indexer copying from input indexer.
 *
 * @param indexer Indexer to copy from.
 */
      RowMajorIndexerBase(const RowMajorIndexerBase<NDIM>& indexer)
      {
        for (unsigned i=0; i<NDIM; ++i)
        {
          start[i] = indexer.start[i];
          shape[i] = indexer.shape[i];
          ai[i] = indexer.ai[i];
        }
        ai[NDIM] = indexer.ai[NDIM];
      }

/**
 * Copy the values from the input indexer.
 *
 * @param indexer Indexer to copy from.
 * @return reference to this indexer.
 */
      RowMajorIndexerBase<NDIM>& operator=(const RowMajorIndexerBase<NDIM>& indexer)
      {
        if (&indexer == this)
          return *this;

        for (unsigned i=0; i<NDIM; ++i)
        {
          start[i] = indexer.start[i];
          shape[i] = indexer.shape[i];
          ai[i] = indexer.ai[i];
        }
        ai[NDIM] = indexer.ai[NDIM];

        return *this;
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
      int getGenIndex(const int idx[NDIM]) const
      {
        int sum = ai[0];
        for (unsigned i=1; i<NDIM+1; ++i)
          sum += ai[i]*idx[i-1];
        return sum;
      }

/**
 * Return linear index given N-dimensional index.
 *
 * @return Linear index.
 */
      int getGenLowIndex(const int idx[NDIM-1]) const
      {
        int sum = ai[0];
        for (unsigned i=1; i<NDIM; ++i)
          sum += ai[i]*idx[i-1];
        sum += ai[NDIM]*start[NDIM-1];
        return sum;
      }

    protected:
/** Coefficients for linear map */
      int ai[NDIM+1];
/** Start indices */
      int start[NDIM];

    private:
/** Shape of linear-space */
      unsigned shape[NDIM];

/**
 * Create a new indexer object.
 *
 * @param shape Shape of region to index.
 * @param start Start index.
 */
      void createRowMajorIndexer(const unsigned shape[NDIM], const int start[NDIM])
      {
        for (unsigned i=0; i<NDIM; ++i)
        {
          this->start[i] = start[i];
          this->shape[i] = shape[i];
        }

        ai[NDIM] = 1;
        for (unsigned i=NDIM-1; i>=1; --i)
          ai[i] = ai[i+1]*shape[i];

        int sum = 0;
        for (unsigned i=1; i<NDIM+1; ++i)
          sum += ai[i]*start[i-1];
        ai[0] = -sum;
      }
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
      RowMajorIndexer(const unsigned shape[NDIM], const int start[NDIM])
        : RowMajorIndexerBase<NDIM>(shape, start)
      {
      }

/**
 * Create a new indexer over given N-dimensional region.
 *
 * @param rgn Region to index.
 */
      RowMajorIndexer(const Region<NDIM, int>& rgn)
        : RowMajorIndexerBase<NDIM>(rgn)
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
      RowMajorIndexer(const unsigned shape[1], const int start[1])
        : RowMajorIndexerBase<1>(shape, start)
      {
      }

/**
 * Create a new indexer over given N-dimensional region.
 *
 * @param rgn Region to index.
 */
      RowMajorIndexer(const Region<1, int>& rgn)
        : RowMajorIndexerBase<1>(rgn)
      {
      }

/**
 * Map 1D index to a linear index.
 *
 * @param i Index location.
 * @return Index of (i) into linear space.
 */
      int getIndex(int i) const 
      {
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
      RowMajorIndexer(const unsigned shape[2], const int start[2])
        : RowMajorIndexerBase<2>(shape, start)
      {
      }

/**
 * Create a new indexer over given N-dimensional region.
 *
 * @param rgn Region to index.
 */
      RowMajorIndexer(const Region<2, int>& rgn)
        : RowMajorIndexerBase<2>(rgn)
      {
      }

/**
 * Map 2D index to a linear index.
 *
 * @param i Index location.
 * @param j Index location.
 * @return Index of (i,j) into linear space.
 */
      int getIndex(int i, int j) const 
      {
        return ai[0]+ai[1]*i+j;
      }

/**
 * Map 2D index to a linear index. The mapping is performed assuming
 * that the final lower indices are the start indices.
 *
 * @param i Index location.
 * @return Index of (i,j) into linear space.
 */
      int getLowIndex(int i) const 
      {
        return getIndex(i, start[1]);
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
      RowMajorIndexer(const unsigned shape[3], const int start[3])
        : RowMajorIndexerBase<3>(shape, start)
      {
      }

/**
 * Create a new indexer over given N-dimensional region.
 *
 * @param rgn Region to index.
 */
      RowMajorIndexer(const Region<3, int>& rgn)
        : RowMajorIndexerBase<3>(rgn)
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
      int getIndex(int i, int j, int k) const 
      {
        return ai[0]+ai[1]*i+ai[2]*j+k;
      }

/**
 * Map 3D index to a linear index. Assume last indices are start
 * indices.
 *
 * @param i Index location.
 * @param j Index location.
 * @return Index of (i,j,s3) into linear space.
 */
      int getLowIndex(int i, int j) const 
      {
        return getIndex(i, j, start[2]);
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
      RowMajorIndexer(const unsigned shape[4], const int start[4])
        : RowMajorIndexerBase<4>(shape, start)
      {
      }

/**
 * Create a new indexer over given N-dimensional region.
 *
 * @param rgn Region to index.
 */
      RowMajorIndexer(const Region<4, int>& rgn)
        : RowMajorIndexerBase<4>(rgn)
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
      int getIndex(int i, int j, int k, int l) const 
      {
        return ai[0]+ai[1]*i+ai[2]*j+ai[3]*k+l;
      }

/**
 * Map 4D index to a linear index.
 *
 * @param i Index location.
 * @param j Index location.
 * @param k Index location.
 * @return Index of (i,j,k,s4) into linear space.
 */
      int getLowIndex(int i, int j, int k) const 
      {
        return getIndex(i, j, k, start[3]);
      }
  };
}

#endif // LC_ROW_MAJOR_INDEXER_H
