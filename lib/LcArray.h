/**
 * @file	LcArray.h
 *
 * @brief	Serial array class.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

#ifndef LC_ARRAY_H
#define LC_ARRAY_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcColMajorIndexer.h>

namespace Lucee
{
// Masks for array traits
  static unsigned int arrayMasks[] =
  {0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80};

// set of macros for setting/getting array traits
#define LC_CONTIGUOUS arrayMasks[0]
#define LC_SET_CONTIGUOUS(bit) (bit) |= LC_CONTIGUOUS
#define LC_CLEAR_CONTIGUOUS(bit) (bit) &= ~LC_CONTIGUOUS
#define LC_IS_CONTIGUOUS(bit) (bit) & LC_CONTIGUOUS

#define LC_ALLOC arrayMasks[1]
#define LC_SET_ALLOC(bit) (bit) |= LC_ALLOC
#define LC_CLEAR_ALLOC(bit) (bit) &= ~LC_ALLOC
#define LC_IS_ALLOC(bit) (bit) & LC_ALLOC

  template <unsigned NDIM, typename T, typename INDEXER = Lucee::ColMajorIndexer<NDIM> >
  class Array
  {
    public:
/**
 * Construct array with specified shape.
 *
 * @param shape Shape of the array.
 * @param init Initial value to assign to all elements.
 */      
      Array(unsigned shape[NDIM], const T& init=(T)0);

/**
 * Construct array with specified shape and start indices.
 *
 * @param shape Shape of the array.
 * @param start indices for the array.
 * @param init Initial value to assign to all elements.
 */
      Array(unsigned shape[NDIM], int start[NDIM], const T& init=(T)0);

/**
 * Destroy object.
 */
      ~Array();

/**
 * Get rank of array.
 *
 * @return Rank of array.
 */
      unsigned getRank() const { return NDIM; }

/**
 * Get total size of array.
 *
 * @return size of array.
 */
      unsigned getSize() const { return len; }

/**
 * Get shape of array.
 *
 * @param shape On return contains the shape of the array.
 */
      void fillWithShape(unsigned shape[NDIM]) const;

/**
 * Get shape of array in the specified direction.
 *
 * @return shape in specified direction.
 */
      unsigned getShape(unsigned dir) const;

/**
 * Get start index in the specified direction.
 *
 * @return start index in specified direction.
 */
      int getLower(unsigned dir) const;

/**
 * Get one past the end index in the specified direction.
 *
 * @return one past the end index in specified direction.
 */
      int getUpper(unsigned dir) const;

/**
 * Is the aray contiguous?
 *
 * @return True if array is contiguous, false otherwise.
 */
      bool isContiguous() const { return LC_IS_CONTIGUOUS(traits); }

/**
 * Accessor function for array.
 *
 * @param idx Indices into array.
 * @return Reference to value at index.
 */
      T& operator()(int idx[NDIM]);

/**
 * Accessor function for array.
 *
 * @param idx Indices into array.
 * @return Value at index.
 */
      T operator()(int idx[NDIM]) const;

/**
 * Accessor function for 1D array.
 *
 * @param i Index into array.
 * @return Reference to value at (i).
 */
      T& operator()(int i);

/**
 * Accessor function for 1D array.
 *
 * @param i Index into array.
 * @return Value at (i).
 */
      T operator()(int i) const;

/**
 * Accessor function for 2D array.
 *
 * @param i Index into array.
 * @param j Index into array.
 * @return Reference to value at (i,j).
 */
      T& operator()(int i, int j);

/**
 * Accessor function for 2D array.
 *
 * @param i Index into array.
 * @param j Index into array.
 * @return Value at (i,j).
 */
      T operator()(int i, int j) const;

/**
 * Accessor function for 3D array.
 *
 * @param i Index into array.
 * @param j Index into array.
 * @param k Index into array.
 * @return Reference to value at (i,j,k).
 */
      T& operator()(int i, int j, int k);

/**
 * Accessor function for 3D array.
 *
 * @param i Index into array.
 * @param j Index into array.
 * @param k Index into array.
 * @return Value at (i,j,k).
 */
      T operator()(int i, int j, int k) const;

/**
 * Accessor function for 4D array.
 *
 * @param i Index into array.
 * @param j Index into array.
 * @param k Index into array.
 * @param l Index into array.
 * @return Reference to value at (i,j,k,l).
 */
      T& operator()(int i, int j, int k, int l);

/**
 * Accessor function for 4D array.
 *
 * @param i Index into array.
 * @param j Index into array.
 * @param k Index into array.
 * @param l Index into array.
 * @return Value at (i,j,k,l).
 */
      T operator()(int i, int j, int k, int l) const;

    private:
/** Indexer into array */
      INDEXER indexer;
/** Array traits stored as bit values */
      unsigned traits;
/** Shape of array */
      unsigned shape[NDIM];
/** Start index of array */
      int start[NDIM];
/** Total size of array */
      unsigned len;
/** Pointer to actual data */
      T *data;
/** Number of arrays pointing to this array */
      mutable int* useCount;

/**
 * Copy constructor is private.
 */
      Array(const Array&);

/**
 * Assignment operator is private.
 */
      Array& operator=(const Array&);

/**
 * This little private class makes a constant array of zeros.
 */
      class Zeros
      {
        public:
          Zeros() {
            for (unsigned i=0; i<NDIM; ++i)
              zeros[i] = 0;
          }
          int zeros[NDIM];
      };
  };

  template <unsigned NDIM, typename T, typename INDEXER>
  Array<NDIM, T, INDEXER>::Array(unsigned shp[NDIM], const T& init)
    : indexer(typename Array<NDIM, T, INDEXER>::Zeros().zeros, shp), traits(0),
      useCount(new int(1))
  {
    len = 1;
    for (unsigned i=0; i<NDIM; ++i)
    {
      shape[i] = shp[i];
      start[i] = 0;
      len = len*shape[i];
    }
    data = new T[len];
    LC_SET_CONTIGUOUS(traits);
    LC_SET_ALLOC(traits);
  }

  template <unsigned NDIM, typename T, typename INDEXER>
  Array<NDIM, T, INDEXER>::Array(unsigned shp[NDIM], int sta[NDIM], const T& init)
    : indexer(sta, shp), traits(0),
      useCount(new int(1))
  {
    len = 1;
    for (unsigned i=0; i<NDIM; ++i)
    {
      shape[i] = shp[i];
      start[i] = sta[i];
      len = len*shape[i];
    }
    data = new T[len];
    LC_SET_CONTIGUOUS(traits);
    LC_SET_ALLOC(traits);
  }

  template <unsigned NDIM, typename T, typename INDEXER>
  Array<NDIM, T, INDEXER>::~Array()
  {
    if ((--*useCount == 0) && LC_IS_ALLOC(traits)) {
      delete [] data;
      delete useCount;
    }
  }

  template <unsigned NDIM, typename T, typename INDEXER>
  void 
  Array<NDIM, T, INDEXER>::fillWithShape(unsigned shp[NDIM]) const
  {
    for (unsigned i=0; i<NDIM; ++i)
      shp[i] = shape[i];
  }

  template <unsigned NDIM, typename T, typename INDEXER>
  unsigned 
  Array<NDIM, T, INDEXER>::getShape(unsigned dir) const
  {
    return shape[dir];
  }

  template <unsigned NDIM, typename T, typename INDEXER>
  int
  Array<NDIM, T, INDEXER>::getLower(unsigned dir) const
  {
    return start[dir];
  }

  template <unsigned NDIM, typename T, typename INDEXER>
  int
  Array<NDIM, T, INDEXER>::getUpper(unsigned dir) const
  {
    return start[dir]+shape[dir];
  }

  template <unsigned NDIM, typename T, typename INDEXER>
  T&
  Array<NDIM, T, INDEXER>::operator()(int idx[NDIM])
  {
    return data[indexer.getGenIndex(idx)];
  }

  template <unsigned NDIM, typename T, typename INDEXER>
  T 
  Array<NDIM, T, INDEXER>::operator()(int idx[NDIM]) const
  {
    return data[indexer.getGenIndex(idx)];
  }

  template <unsigned NDIM, typename T, typename INDEXER>
  T&
  Array<NDIM, T, INDEXER>::operator()(int i)
  {
    return data[indexer.getIndex(i)];
  }

  template <unsigned NDIM, typename T, typename INDEXER>
  T
  Array<NDIM, T, INDEXER>::operator()(int i) const
  {
    return data[indexer.getIndex(i)];
  }

  template <unsigned NDIM, typename T, typename INDEXER>
  T&
  Array<NDIM, T, INDEXER>::operator()(int i, int j)
  {
    return data[indexer.getIndex(i, j)];
  }

  template <unsigned NDIM, typename T, typename INDEXER>
  T
  Array<NDIM, T, INDEXER>::operator()(int i, int j) const
  {
    return data[indexer.getIndex(i, j)];
  }

  template <unsigned NDIM, typename T, typename INDEXER>
  T&
  Array<NDIM, T, INDEXER>::operator()(int i, int j, int k)
  {
    return data[indexer.getIndex(i, j, k)];
  }

  template <unsigned NDIM, typename T, typename INDEXER>
  T
  Array<NDIM, T, INDEXER>::operator()(int i, int j, int k) const
  {
    return data[indexer.getIndex(i, j, k)];
  }

  template <unsigned NDIM, typename T, typename INDEXER>
  T& 
  Array<NDIM, T, INDEXER>::operator()(int i, int j, int k, int l)
  {
    return data[indexer.getIndex(i, j, k, l)];
  }

  template <unsigned NDIM, typename T, typename INDEXER>
  T 
  Array<NDIM, T, INDEXER>::operator()(int i, int j, int k, int l) const
  {
    return data[indexer.getIndex(i, j, k, l)];
  }
}

#endif // LC_ARRAY_H
