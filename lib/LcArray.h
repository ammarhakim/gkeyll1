/**
 * @file	LcArray.h
 *
 * @brief	Serial array class.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_ARRAY_H
#define LC_ARRAY_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcColMajorIndexer.h>
#include <LcDataStructIfc.h>
#include <LcExcept.h>
#include <LcFixedVector.h>
#include <LcRegion.h>
#include <LcRowMajorSequencer.h>

// std includes
#include <vector>
#include <string>

namespace Lucee
{
// Masks for array traits
  static unsigned int arrayMasks[] =
  {0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80};

// set of macros for setting/getting array traits

/** Contiguous mask */
#define LC_CONTIGUOUS arrayMasks[0]
/** Set contiguous mask */
#define LC_SET_CONTIGUOUS(bit) (bit) |= LC_CONTIGUOUS
/** Clear contiguous mask */
#define LC_CLEAR_CONTIGUOUS(bit) (bit) &= ~LC_CONTIGUOUS
/** Check if contiguous mask is set */
#define LC_IS_CONTIGUOUS(bit) (bit) & LC_CONTIGUOUS

/** Allocation mask */
#define LC_ALLOC arrayMasks[1]
/** Set allocation mask */
#define LC_SET_ALLOC(bit) (bit) |= LC_ALLOC
/** Clear allocation mask */
#define LC_CLEAR_ALLOC(bit) (bit) &= ~LC_ALLOC
/** Check if allocation mask is set */
#define LC_IS_ALLOC(bit) (bit) & LC_ALLOC

/**
 * Basic multi-dimensional array class.
 */
  template <unsigned NDIM, typename T, template <unsigned> class INDEXER  = Lucee::ColMajorIndexer>
  class Array : public DataStructIfc
  {
    public:
/** We need to friend ourself to allow accessing private stuff from another dimension */
      template <unsigned RDIM, typename TT, template <unsigned> class IINDEXER> friend class Array;

/**
 * Construct array with specified shape.
 *
 * @param shape Shape of the array.
 * @param init Initial value to assign to all elements.
 */
      Array(const unsigned shape[NDIM], const T& init=(T)0);

/**
 * Construct array with specified shape and start indices.
 *
 * @param shape Shape of the array.
 * @param start Start indices for the array.
 * @param init Initial value to assign to all elements.
 */
      Array(const unsigned shape[NDIM], const int start[NDIM], const T& init=(T)0);

/**
 * Construct array with specified index region.
 *
 * @param rgn Region indexed by array.
 * @param init Initial value to assign to all elements.
 */
      Array(const Lucee::Region<NDIM, int>& rgn, const T& init=(T)0);

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
 * Make a copy of the array.
 *
 * @return copy of the array.
 */
      Lucee::Array<NDIM, T> duplicate();

/**
 * Return reference to the first element in array.
 *
 * @return Reference to first element in array.
 */
      T& first();

/**
 * Region indexed by array.
 *
 * @return region indexed by array.
 */
      Lucee::Region<NDIM, int> getRegion() const;

/**
 * Return a slice into the array. The returned array shares the data
 * with this array.
 *
 * @param rgn Region of the slice.
 * @return Sliced array.
 */
      Array<NDIM, T, INDEXER> getSlice(const Lucee::Region<NDIM, int>& rgn);

/**
 * Deflate the array into a lower-dimensional array. The returned
 * array shares the data with this array.
 *
 * @param defDims Dimensions to remove from array.
 * @param defDimsIdx Index along the removed dimensions.
 * @return deflated array.
 */
      template <unsigned RDIM>
      Array<RDIM, T, INDEXER> 
      deflate(const unsigned defDims[NDIM-RDIM], const int defDimsIdx[NDIM-RDIM])
      {
        ++*useCount; // increment use count in case array is deleted before slice
// extend defDims array for use in creating new lower/upper bounds
        unsigned extDefDims[NDIM];
        for (unsigned i=0; i<NDIM-RDIM; ++i)
          extDefDims[i] = defDims[i];
        for (unsigned i=NDIM-RDIM; i<NDIM; ++i)
          extDefDims[i] = NDIM+1; // so that cmp in if-statement fails

// construct deflated lower and upper bounds
        int defLower[RDIM], defUpper[RDIM];
        unsigned ddIdx = 0, nddIdx = 0;
        for (unsigned i=0; i<NDIM; ++i)
        {
          if (extDefDims[ddIdx] != i)
          {
            defLower[nddIdx] = getLower(i);
            defUpper[nddIdx] = getUpper(i);
            nddIdx++;
          }
          else
          {
            ddIdx++;
          }
        }
        Lucee::Region<RDIM, int> defRgn(defLower, defUpper);
        Lucee::Array<RDIM, T, INDEXER> na(defRgn, data);
        na.indexer = indexer.template deflate<RDIM>(defDims, defDimsIdx);
        delete na.useCount;
        na.useCount = useCount;
        
        return na;
      }

/**
 * Accessor function for array.
 *
 * @param idx Indices into array.
 * @return Reference to value at index.
 */
      T& operator()(const int idx[NDIM]);

/**
 * Accessor function for array.
 *
 * @param idx Indices into array.
 * @return Value at index.
 */
      T operator()(const int idx[NDIM]) const;

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

/**
 * Reset the lower bounds of the array.
 *
 * @param nlo New lower bounds.
 */
      void resetLower(const int nlo[NDIM]);

/**
 * Write dataStruct to given node in HDF5 file.
 *
 * @param io I/O object for I/O.
 * @param node Node to write to.
 * @param nm Name of the grid as it should appear in output.
 * @return node to which data was written.
 */
      virtual Lucee::IoNodeType writeToFile(Lucee::IoBase& io, Lucee::IoNodeType& node,
        const std::string& nm);

    protected:
/**
 * Construct array with specified index region and re-using the
 * supplied data-space.
 *
 * @param rgn Region indexed by array.
 * @param dp Pointer to data space to use.
 */
      Array(const Lucee::Region<NDIM, int>& rgn, T* dp);

/**
 * Returns reference to location in underlying memory space.
 *
 * @param loc Location into memory space.
 * @return reference to data.
 */
      T& getRefToLoc(unsigned loc);

/**
 * Returns reference to location in underlying memory space.
 *
 * @param loc Location into memory space.
 * @return reference to data.
 */
      const T& getConstRefToLoc(unsigned loc) const;

    private:
/** Indexer into array */
      INDEXER<NDIM> indexer;
/** Array traits stored as bit values */
      unsigned traits;
/** Shape of array */
      int shape[NDIM];
/** Start index of array */
      int start[NDIM];
/** Total size of array */
      unsigned len;
/** Number of arrays pointing to this array */
      mutable int* useCount;
/** Pointer to actual data */
      T *data;
  };
  
  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  Array<NDIM, T, INDEXER>::Array(const unsigned shp[NDIM], const T& init)
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
    for (unsigned i=0; i<len; ++i)
      data[i] = init;

    LC_SET_CONTIGUOUS(traits);
    LC_SET_ALLOC(traits);
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  Array<NDIM, T, INDEXER>::Array(const unsigned shp[NDIM], const int sta[NDIM], const T& init)
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
    for (unsigned i=0; i<len; ++i)
      data[i] = init;

    LC_SET_CONTIGUOUS(traits);
    LC_SET_ALLOC(traits);
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  Array<NDIM, T, INDEXER>::Array(const Lucee::Region<NDIM, int>& rgn, const T& init)
    : indexer(rgn), traits(0), useCount(new int(1))
  {
    len = 1;
    for (unsigned i=0; i<NDIM; ++i)
    {
      shape[i] = rgn.getShape(i);
      start[i] = rgn.getLower(i);
      len = len*shape[i];
    }
    data = new T[len];
    for (unsigned i=0; i<len; ++i)
      data[i] = init;

    LC_SET_CONTIGUOUS(traits);
    LC_SET_ALLOC(traits);
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
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

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
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

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  Array<NDIM, T, INDEXER>::~Array()
  {
    if ((--*useCount == 0) && LC_IS_ALLOC(traits)) 
    {
      delete [] data;
      delete useCount;
    }
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  Array<NDIM, T, INDEXER>&
  Array<NDIM, T, INDEXER>::operator=(const T& val)
  {
    if (isContiguous())
// just set value directly into data-space
      for (unsigned i=0; i<getSize(); ++i)
        data[i] = val;
    else
    {
// create sequencer
      Lucee::Region<NDIM, int> rgn
        = Lucee::createRegionFromStartAndShape<NDIM, int>(start, shape);
      typename INDEXER<NDIM>::Sequencer seq(rgn);
// loop over region
      while (seq.step())
        data[indexer.getIndex(seq.getIndex())] = val;
    }

    return *this;
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  void 
  Array<NDIM, T, INDEXER>::fillWithShape(unsigned shp[NDIM]) const
  {
    for (unsigned i=0; i<NDIM; ++i)
      shp[i] = shape[i];
  }


  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  void 
  Array<NDIM, T, INDEXER>::fillWithStart(int strt[NDIM]) const
  {
    for (unsigned i=0; i<NDIM; ++i)
      strt[i] = start[i];
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  unsigned 
  Array<NDIM, T, INDEXER>::getShape(unsigned dir) const
  {
    return shape[dir];
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  int
  Array<NDIM, T, INDEXER>::getLower(unsigned dir) const
  {
    return start[dir];
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  int
  Array<NDIM, T, INDEXER>::getUpper(unsigned dir) const
  {
    return start[dir]+shape[dir];
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  Lucee::Array<NDIM, T>
  Array<NDIM, T, INDEXER>::duplicate()
  {
    Lucee::Region<NDIM, int> rgn
      = Lucee::createRegionFromStartAndShape<NDIM, int>(start, shape);
    Lucee::Array<NDIM, T, INDEXER> dup(rgn);
// create sequencer to copy stuff over
    typename INDEXER<NDIM>::Sequencer seq(rgn);
// copy into new array
    while (seq.step())
      dup(seq.getIndex()) = this->operator()(seq.getIndex());

    return dup;
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  T&
  Array<NDIM, T, INDEXER>::first()
  {
    return data[indexer.getIndex(start)];
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  Lucee::Region<NDIM, int>
  Array<NDIM, T, INDEXER>::getRegion() const
  {
    int upper[NDIM];
    for (unsigned i=0; i<NDIM; ++i)
      upper[i] = start[i]+shape[i];
    return Lucee::Region<NDIM, int>(start, upper);
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  Array<NDIM, T, INDEXER>
  Array<NDIM, T, INDEXER>::getSlice(const Lucee::Region<NDIM, int>& rgn)
  {
// check if specified region is contained in this    
    if (!getRegion().contains(rgn))
      throw Lucee::Except("Array::getSlice: Slice region must be contained in the array region");

    ++*useCount; // increment use count in case array is deleted before slice
    Lucee::Array<NDIM, T, INDEXER> na(rgn, data);
    na.indexer = indexer; // slice uses same indexer as this array
    delete na.useCount;
    na.useCount = useCount;

    return na;
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  T&
  Array<NDIM, T, INDEXER>::operator()(const int idx[NDIM])
  {
    return data[indexer.getIndex(idx)];
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  T 
  Array<NDIM, T, INDEXER>::operator()(const int idx[NDIM]) const
  {
    return data[indexer.getIndex(idx)];
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  T&
  Array<NDIM, T, INDEXER>::operator()(int i)
  {
    return data[indexer.getIndex(i)];
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  T
  Array<NDIM, T, INDEXER>::operator()(int i) const
  {
    return data[indexer.getIndex(i)];
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  T&
  Array<NDIM, T, INDEXER>::operator()(int i, int j)
  {
    return data[indexer.getIndex(i, j)];
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  T
  Array<NDIM, T, INDEXER>::operator()(int i, int j) const
  {
    return data[indexer.getIndex(i, j)];
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  T&
  Array<NDIM, T, INDEXER>::operator()(int i, int j, int k)
  {
    return data[indexer.getIndex(i, j, k)];
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  T
  Array<NDIM, T, INDEXER>::operator()(int i, int j, int k) const
  {
    return data[indexer.getIndex(i, j, k)];
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  T& 
  Array<NDIM, T, INDEXER>::operator()(int i, int j, int k, int l)
  {
    return data[indexer.getIndex(i, j, k, l)];
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  T 
  Array<NDIM, T, INDEXER>::operator()(int i, int j, int k, int l) const
  {
    return data[indexer.getIndex(i, j, k, l)];
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  T& 
  Array<NDIM, T, INDEXER>::getRefToLoc(unsigned loc)
  {
    return data[loc];
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  const T& 
  Array<NDIM, T, INDEXER>::getConstRefToLoc(unsigned loc) const
  {
    return data[loc];
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  void
  Array<NDIM, T, INDEXER>::resetLower(const int nlo[NDIM])
  {
    for (unsigned i=0; i<NDIM; ++i)
      start[i] = nlo[i];
    indexer.resetLower(nlo);
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  Array<NDIM, T, INDEXER>::Array(const Lucee::Region<NDIM, int>& rgn, T* dp)
    : indexer(rgn), traits(0), useCount(new int(1)), data(dp)
  {
    len = 1;
    for (unsigned i=0; i<NDIM; ++i)
    {
      shape[i] = rgn.getShape(i);
      start[i] = rgn.getLower(i);
      len = len*shape[i];
    }
    LC_CLEAR_CONTIGUOUS(traits);
    LC_CLEAR_ALLOC(traits);
  }

  template <unsigned NDIM, typename T, template <unsigned> class INDEXER>
  Lucee::IoNodeType
  Array<NDIM, T, INDEXER>::writeToFile(Lucee::IoBase& io, 
    Lucee::IoNodeType& node, const std::string& nm)
  {
    std::vector<size_t> dataSetSize(NDIM), dataSetBeg(NDIM), dataSetLen(NDIM);
// construct sizes and shapes to write stuff out
    for (unsigned i=0; i<NDIM; ++i)
    {
      dataSetSize[i] = this->getShape(i);
      dataSetBeg[i] = 0;
      dataSetLen[i] = this->getShape(i);
    }

    Lucee::Region<NDIM, int> rgn = this->getRegion();
    std::vector<T> buff(rgn.getVolume());
    Lucee::RowMajorSequencer<NDIM> seq(rgn); // must be row-major for HDF5
// copy data into buffer
    unsigned count = 0;
    while (seq.step())
      buff[count++] = this->operator()(seq.getIndex());
// write it out
    Lucee::IoNodeType dn =
      io.writeDataSet(node, nm, dataSetSize, dataSetBeg, dataSetLen, &buff[0]);

    return dn;
  }
}

#endif // LC_ARRAY_H
