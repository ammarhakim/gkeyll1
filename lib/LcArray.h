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
#include <LcExcept.h>
#include <LcFixedVector.h>

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
 * @param start Start indices for the array.
 * @param init Initial value to assign to all elements.
 */
      Array(unsigned shape[NDIM], int start[NDIM], const T& init=(T)0);

/**
 * Create this array from the input array. The copy does not allocate
 * new space but simply points to the input array.
 *
 * @param arr Array to create from.
 */
      Array(const Array<NDIM, T, INDEXER>& arr);

/**
 * Copy array from input array. The copy does not allocate new space
 * but simply points to the input array.
 *
 * @param arr Array to copy from.
 * @return Reference to this array.
 */
      Array<NDIM, T, INDEXER>& operator=(const Array<NDIM, T, INDEXER>& arr);

/**
 * Destroy object.
 */
      ~Array();

/**
 * Assign all elements in array to specified value.
 *
 * @param val Value to assign.
 * @return Reference to this array.
 */
      Array<NDIM, T, INDEXER>& operator=(const T& val);

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
 * Get start index of array.
 *
 * @param start On return contains the start index of array.
 */
      void fillWithStart(int start[NDIM]) const;

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
 * @return One past the end index in specified direction.
 */
      int getUpper(unsigned dir) const;

/**
 * Is the array contiguous?
 *
 * @return True if array is contiguous, false otherwise.
 */
      bool isContiguous() const { return LC_IS_CONTIGUOUS(traits); }

/**
 * Return reference to the first element in array.
 *
 * @param Reference to first element in array.
 */
      T& first();

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

    protected:
/**
 * Create a new array without specifying its shape or start
 * indices. Data is allocated but region is not specified.
 *
 * @param nelem Number of elemenets in array.
 */
      Array(unsigned nelem);

/**
 * Create a new array without specifying its shape or start
 * indices. Data is not allocated, the supplied pointer is used
 * instead. No checks are done to ensure that supplied pointer has
 * sufficient space.
 *
 * @param nelem Number of elemenets in array.
 * @param dataSpace Data space to re-use.
 */
      Array(unsigned nelem, T *dataSpace);

/**
 * Set array shape and start indices.
 *
 * @param shape Shape of array.
 * @param start Start indices of array.
 */
      void setRegion(unsigned shape[NDIM], int start[NDIM]);

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
  };

  template <unsigned NDIM, typename T, typename INDEXER>
  Array<NDIM, T, INDEXER>::Array(unsigned shp[NDIM], const T& init)
    : indexer(shp, &Lucee::FixedVector<NDIM,int>(0)[0]),
      traits(0),
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
    : indexer(shp, sta), traits(0),
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
  Array<NDIM, T, INDEXER>::Array(const Array<NDIM, T, INDEXER>& arr)
    : indexer(arr.indexer), traits(arr.traits), len(arr.len), data(arr.data)
  {
    for (unsigned i=0; i<NDIM; ++i)
    {
      start[i] = arr.start[i];
      shape[i] = arr.shape[i];
    }
    LC_CLEAR_ALLOC(traits); // no memory was allocated
// increment use-count of input array
    ++*arr.useCount;
    useCount = arr.useCount;
  }

  template <unsigned NDIM, typename T, typename INDEXER>
  Array<NDIM, T, INDEXER>&
  Array<NDIM, T, INDEXER>::operator=(const Array<NDIM, T, INDEXER>& arr)
  {
    if (&arr == this)
      return *this;

// increment use count of input array
    ++*arr.useCount; 

// decrement our use count and delete memory if no one else is
// pointing to us
    if (--*useCount == 0)
    {
      if (LC_IS_ALLOC(traits))
      {
        delete [] data;
        delete useCount;
      }
    }

// copy stuff over
    indexer = arr.indexer;
    traits = arr.traits;
    len = arr.len;

    for (unsigned i=0; i<NDIM; ++i)
    {
      start[i] = arr.start[i];
      shape[i] = arr.shape[i];
    }
    LC_CLEAR_ALLOC(traits); // no memory was allocated

    data = arr.data;
    useCount = arr.useCount;

    return *this;
  }

  template <unsigned NDIM, typename T, typename INDEXER>
  Array<NDIM, T, INDEXER>::~Array()
  {
    if ((--*useCount == 0) && LC_IS_ALLOC(traits)) 
    {
      delete [] data;
      delete useCount;
    }
  }

  template <unsigned NDIM, typename T, typename INDEXER>
  Array<NDIM, T, INDEXER>&
  Array<NDIM, T, INDEXER>::operator=(const T& val)
  {
    if (isContiguous())
// just copy value directly into data-space
      for (unsigned i=0; i<getSize(); ++i)
        data[i] = val;
    else
    {
      throw Lucee::Except("Array::operator=: NOT IMPLEMENTED");
    }

    return *this;
  }

  template <unsigned NDIM, typename T, typename INDEXER>
  void 
  Array<NDIM, T, INDEXER>::fillWithShape(unsigned shp[NDIM]) const
  {
    for (unsigned i=0; i<NDIM; ++i)
      shp[i] = shape[i];
  }


  template <unsigned NDIM, typename T, typename INDEXER>
  void 
  Array<NDIM, T, INDEXER>::fillWithStart(int strt[NDIM]) const
  {
    for (unsigned i=0; i<NDIM; ++i)
      strt[i] = start[i];
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
  Array<NDIM, T, INDEXER>::first()
  {
    if (isContiguous() == false)
      throw Lucee::Except("Array::first: Array must be contiguous.");
    return data[indexer.getGenIndex(start)];
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

  template <unsigned NDIM, typename T, typename INDEXER>
  Array<NDIM, T, INDEXER>::Array(unsigned nelem)
    : len(nelem), data(new T[nelem]), useCount(new int(1))
  {
    LC_SET_CONTIGUOUS(traits);
    LC_SET_ALLOC(traits);
  }

  template <unsigned NDIM, typename T, typename INDEXER>
  Array<NDIM, T, INDEXER>::Array(unsigned nelem, T *dataSpace)
    : len(nelem), data(dataSpace), useCount(new int(1))
  {
    LC_SET_CONTIGUOUS(traits);
    LC_CLEAR_ALLOC(traits); // we did not allocate this array
  }

  template <unsigned NDIM, typename T, typename INDEXER>
  void
  Array<NDIM, T, INDEXER>::setRegion(unsigned shp[NDIM], int sta[NDIM])
  {
    for (unsigned i=0; i<NDIM; ++i)
    {
      shape[i] = shp[i];
      start[i] = sta[i];
    }
  }
}

#endif // LC_ARRAY_H
