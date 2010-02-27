/**
 * @file	LcLinIndexer.h
 *
 * @brief	Base class mapping an N-dimensional index into a linear index.
 *
 * @version	$Id: LcLinIndexer.h 280 2010-02-22 23:44:44Z a.hakim777 $
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

#ifndef LC_LIN_INDEXER_H
#define LC_LIN_INDEXER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcFixedVector.h>
#include <LcRegion.h>

namespace Lucee
{
/**
 * Base class provides key coefficient contruction functionality for
 * use in the indexing functions provided by the derived classes.
 */
  template <unsigned NDIM>
  class LinIndexer
  {
    public:
/**
 * Create a new indexer for mapping an N-dimensional index into a
 * linear index using a linear mapping.
 *
 * @param shape Shape of space.
 * @param start Starting index.
 * @param coeff Coefficients for linear mapping.
 */
      LinIndexer(const Lucee::FixedVector<NDIM, unsigned>& shape,
        const Lucee::FixedVector<NDIM, int>& start, const Lucee::FixedVector<NDIM+1, int>& coeff);

/**
 * Create a new indexer for mapping an N-dimensional index into a
 * linear index using a linear mapping.
 *
 * @param rgn Region for space.
 * @param coeff Coefficients for linear mapping.
 */
      LinIndexer(const Lucee::Region<NDIM, int>& rgn, const Lucee::FixedVector<NDIM+1, int>& coeff);

/**
 * Create a new indexer copying from input indexer.
 *
 * @param indexer Indexer to copy from.
 */
      LinIndexer(const LinIndexer<NDIM>& indexer);

/**
 * Copy the values from the input indexer.
 *
 * @param indexer Indexer to copy from.
 * @return reference to this indexer.
 */
      LinIndexer<NDIM>& operator=(const LinIndexer<NDIM>& indexer);

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
 * Map 1D index to a linear index.
 *
 * @param i Index location.
 * @return Index of (i) into linear space.
 */
      int getIndex(int i) const;

/**
 * Map 2D index to a linear index.
 *
 * @param i Index location.
 * @param j Index location.
 * @return Index of (i, j) into linear space.
 */
      int getIndex(int i, int j) const;

/**
 * Map 3D index to a linear index.
 *
 * @param i Index location.
 * @param j Index location.
 * @param k Index location.
 * @return Index of (i, j, k) into linear space.
 */
      int getIndex(int i, int j, int k) const;

/**
 * Map 4D index to a linear index.
 *
 * @param i Index location.
 * @param j Index location.
 * @param k Index location.
 * @param l Index location.
 * @return Index of (i, j, k, l) into linear space.
 */
      int getIndex(int i, int j, int k, int l) const;

/**
 * Map NDIM index to a linear index.
 *
 * @param idx NDIM-D index.
 * @return Index into linear space.
 */
      int getGenIndex(const int idx[NDIM]) const;

/**
 * Map 2D index to a linear index.
 *
 * @param i Index location.
 * @return Index of (i, s2) into linear space.
 */
      int getLowIndex(int i) const;

/**
 * Map 3D index to a linear index.
 *
 * @param i Index location.
 * @param j Index location.
 * @return Index of (i, j, s3) into linear space.
 */
      int getLowIndex(int i, int j) const;

/**
 * Map 4D index to a linear index.
 *
 * @param i Index location.
 * @param j Index location.
 * @param k Index location.
 * @return Index of (i, j, k, s4) into linear space.
 */
      int getLowIndex(int i, int j, int k) const;

/**
 * Map NDIM index to a linear index.
 *
 * @param idx NDIM-D index.
 * @return Index into linear space.
 */
      int getGenLowIndex(const int idx[NDIM-1]) const;

    protected:
/**
 * Return an indexer that allows indexing a lower-dimensional sub-space.
 *
 * @param defDims Dimensions to remove from space.
 * @param defDimsIdx Index along the removed dimensions.
 * @return deflated indexer.
 */
      template <unsigned RDIM>
      LinIndexer<RDIM>
      deflateLin(const unsigned defDims[NDIM-RDIM], const int defDimsIdx[NDIM-RDIM]) const
      {
// extend defDims array for use in creating new lower/upper bounds
        unsigned extDefDims[NDIM];
        for (unsigned i=0; i<NDIM-RDIM; ++i)
          extDefDims[i] = defDims[i];
        for (unsigned i=NDIM-RDIM; i<NDIM; ++i)
          extDefDims[i] = NDIM+1; // so that cmp in if-statement fails

// construct deflated-indexer start, shape and coefficients
        int defStart[RDIM];
        unsigned defShape[RDIM], ddIdx = 0, nddIdx = 0;
        int defAi[RDIM+1];
// compute 0th coefficient
        defAi[0] = ai[0];
        for (unsigned i=0; i<NDIM-RDIM; ++i)
          defAi[0] += defDimsIdx[i]*ai[defDims[i]+1];
// now compute the other stuff
        for (unsigned i=0; i<NDIM; ++i)
        {
          if (extDefDims[ddIdx] != i)
          {
            defStart[nddIdx] = start[i];
            defShape[nddIdx] = shape[i];
            defAi[nddIdx+1] = ai[i+1];
            nddIdx++;
          }
          else
          {
            ddIdx++;
          }
        }
        return LinIndexer<RDIM>(defShape, defStart, defAi);
      }

    private:
/** Start indices */
      Lucee::FixedVector<NDIM, int> start;
/** Shape of linear-space */
      Lucee::FixedVector<NDIM, unsigned> shape;
/** Coefficients for linear map */
      int ai[NDIM+1];
  };

  template <unsigned NDIM>
  LinIndexer<NDIM>::LinIndexer(const Lucee::FixedVector<NDIM, unsigned>& shape,
    const Lucee::FixedVector<NDIM, int>& start, const Lucee::FixedVector<NDIM+1, int>& coeff)
    : start(start), shape(shape)
  {
    for (unsigned i=0; i<NDIM; ++i)
      ai[i] = coeff[i];
    ai[NDIM] = coeff[NDIM];
  }

  template <unsigned NDIM>
  LinIndexer<NDIM>::LinIndexer(const Lucee::Region<NDIM, int>& rgn, const Lucee::FixedVector<NDIM+1, int>& coeff)
    : start(rgn.getLower()), shape((unsigned) 0)
  {
    for (unsigned i=0; i<NDIM; ++i)
    {
      shape[i] = rgn.getShape(i);
      ai[i] = coeff[i];
    }
    ai[NDIM] = coeff[NDIM];
  }

  template <unsigned NDIM>
  LinIndexer<NDIM>::LinIndexer(const LinIndexer<NDIM>& indexer)
    : start(indexer.start), shape(indexer.shape)
  {
    for (unsigned i=0; i<NDIM; ++i)
      ai[i] = indexer.ai[i];
    ai[NDIM] = indexer.ai[NDIM];
  }

  template <unsigned NDIM>
  LinIndexer<NDIM>&
  LinIndexer<NDIM>::operator=(const LinIndexer<NDIM>& indexer)
  {
    if (&indexer == this)
      return *this;

    start = indexer.start;
    shape = indexer.shape;    
    for (unsigned i=0; i<NDIM; ++i)
      ai[i] = indexer.ai[i];
    ai[NDIM] = indexer.ai[NDIM];
    
    return *this;
  }

  template <unsigned NDIM>
  int
  LinIndexer<NDIM>::getIndex(int i) const
  {
    return ai[0] + i*ai[1];
  }

  template <unsigned NDIM>
  int
  LinIndexer<NDIM>::getIndex(int i, int j) const
  {
    return ai[0] + i*ai[1] + j*ai[2];
  }

  template <unsigned NDIM>
  int
  LinIndexer<NDIM>::getIndex(int i, int j, int k) const
  {
    return ai[0] + i*ai[1] + j*ai[2] + k*ai[3];
  }

  template <unsigned NDIM>
  int
  LinIndexer<NDIM>::getIndex(int i, int j, int k, int l) const
  {
    return ai[0] + i*ai[1] + j*ai[2] + k*ai[3] + l*ai[4];
  }

  template <unsigned NDIM>
  int
  LinIndexer<NDIM>::getGenIndex(const int idx[NDIM]) const
  {
    int sum = ai[0];
    for (unsigned i=1; i<NDIM+1; ++i)
      sum += idx[i-1]*ai[i];
    return sum;
  }

  template <unsigned NDIM>
  int
  LinIndexer<NDIM>::getLowIndex(int i) const
  {
    return getIndex(i, start[1]);
  }

  template <unsigned NDIM>
  int
  LinIndexer<NDIM>::getLowIndex(int i, int j) const
  {
    return getIndex(i, j, start[2]);
  }

  template <unsigned NDIM>
  int
  LinIndexer<NDIM>::getLowIndex(int i, int j, int k) const
  {
    return getIndex(i, j, k, start[3]);
  }

  template <unsigned NDIM>
  int
  LinIndexer<NDIM>::getGenLowIndex(const int idx[NDIM-1]) const
  {
    int sum = ai[0];
    for (unsigned i=1; i<NDIM; ++i)
      sum += ai[i]*idx[i-1];
    sum += ai[NDIM]*start[NDIM-1];
    return sum;
  }
}

#endif // LC_LIN_INDEXER_H
