/**
 * @file	LcField.h
 *
 * @brief	Fields hold multiple values per index location.
 *
 * @version	$Id: LcField.h 222 2009-11-17 04:46:22Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_FIELD_H
#define LC_FIELD_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcArray.h>
#include <LcConstFieldPtr.h>
#include <LcFieldPtr.h>
#include <LcRegion.h>
#include <LcRowMajorIndexer.h>

namespace Lucee
{
/**
 * A field represents an array that can hold multiple values per index
 * location. Fields can be indexed directly using (i,j,k,..,c)
 * notation, with the last index 'c' is into the field components, or
 * using field iterators. Further, each field can have an extra set of
 * indices outside the main indexed region which can be used to
 * represent "ghost cells or ghost nodes". This is useful, for
 * example, in computing finite-difference stencils on the field.
 */
  template <unsigned NDIM, typename T>
  class Field : public Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer>
  {
    public:
/** We need to friend ourself to allow accessing private stuff from another dimension */
      template <unsigned RDIM, typename TT> friend class Field;

/**
 * Create a new field indexing given region. This constructor creates
 * an empty set of ghost indices.
 *
 * @param rgn Region indexed by array.
 * @param nc Number of components at each index location.
 * @param init Inital value to assigned to all components.
 */      
      Field(const Lucee::Region<NDIM, int>& rgn, unsigned nc, const T& init=(T)0);

/**
 * Create a new field indexing given region. This constructor creates
 * an empty set of ghost indices.
 *
 * @param rgn Region indexed by array.
 * @param nc Number of components at each index location.
 * @param lg Ghost indexes along lower index range in each dimension.
 * @param ug Ghost indexes along upper index range in each dimension.
 * @param init Inital value to assigned to all components.
 */      
      Field(const Lucee::Region<NDIM, int>& rgn, unsigned nc, 
        int lg[NDIM], int ug[NDIM], const T& init=(T)0);

/**
 * Create a field from supplied one (shallow copy).
 *
 * @param fld Field to copy from.
 */
      Field(const Field<NDIM, T>& fld);

/**
 * Assign a field from supplied one (shallow copy).
 *
 * @param fld Field to assign from.
 * @return Reference to this field.
 */
      Field<NDIM, T>& operator=(const Field<NDIM, T>& fld);

/**
 * Assign field to given value.
 *
 * @param val Value to assign.
 * @return reference to this field.
 */
      Field<NDIM, T> & operator=(const T& val);

/**
 * Number of components per index location.
 *
 * @return number of components.
 */
      unsigned getNumComponents() const { return numComponents; }

/**
 * Get start index in the specified direction.
 *
 * @return start index in specified direction.
 */
      int getLower(unsigned dir) const
      { return rgn.inflate(0, numComponents).getLower(dir); }

/**
 * Get one past the end index in the specified direction.
 *
 * @return One past the end index in specified direction.
 */
      int getUpper(unsigned dir) const
      { return rgn.inflate(0, numComponents).getUpper(dir); }

/**
 * Get start index of the extended region (including ghost indices) in
 * the specified direction.
 *
 * @return start index in specified direction.
 */
      int getLowerExt(unsigned dir) const
      { 
        return rgn.extend(lowerGhost, upperGhost)
          .inflate(0, numComponents).getLower(dir); 
      }

/**
 * Get one past the end index of the extended region (including ghost
 * indices) in the specified direction.
 *
 * @return One past the end index in specified direction.
 */
      int getUpperExt(unsigned dir) const
      { 
        return rgn.extend(lowerGhost, upperGhost)
          .inflate(0, numComponents).getUpper(dir);
      }

/**
 * Region indexed by field.
 *
 * @return region indexed by field.
 */
      Lucee::Region<NDIM, int> getRegion() const 
      { return rgn; }

/**
 * Region extended region (including ghost indices) indexed by field.
 *
 * @return extended region indexed by field.
 */
      Lucee::Region<NDIM, int> getExtRegion() const 
      { return rgn.extend(lowerGhost, upperGhost); }

/**
 * Get indexer into field.
 *
 * @return Indexer into field.
 */
      Lucee::RowMajorIndexer<NDIM+1> getIndexer() const { return rgnIdx; }

/**
 * Get a view into the field. The returned object has the
 * same-dimensionality and shares data with this object. The view has
 * access to all components stored in this field.
 *
 * @param rgn Region of the slice.
 * @return View into field.
 */
      Field<NDIM, T> getView(const Lucee::Region<NDIM, int>& rgn);

/**
 * Get a field with same indexed region as this one, but with only a
 * range of the components indexed. The returned field has access to
 * components [sc, ec).
 *
 * @param sc Index into components.
 * @param ec One past the last index into components.
 * @return View into field.
 */
      Field<NDIM, T> getSubCompView(unsigned sc, unsigned ec);

/**
 * Create a new pointer object to elements in field.
 *
 * @return Pointer to first-element in field.
 */
      Lucee::FieldPtr<T> createPtr();

/**
 * Create a new pointer object to elements in field.
 *
 * @return Pointer to first-element in field.
 */
      Lucee::ConstFieldPtr<T> createConstPtr() const;

/**
 * Set pointer to given (i) 1D location.
 *
 * @param ptr Pointer to set.
 * @param i Location to set to.
 */
      void setPtr(Lucee::FieldPtr<T>& ptr, int i);

/**
 * Set pointer to given (i,j) 2D location.
 *
 * @param ptr Pointer to set.
 * @param i Location to set to.
 * @param j Location to set to.
 */
      void setPtr(Lucee::FieldPtr<T>& ptr, int i, int j);

/**
 * Set pointer to given (i,j,k) 3D location.
 *
 * @param ptr Pointer to set.
 * @param i Location to set to.
 * @param j Location to set to.
 * @param k Location to set to.
 */
      void setPtr(Lucee::FieldPtr<T>& ptr, int i, int j, int k);

/**
 * Set pointer to given N-dimensional location.
 *
 * @param ptr Pointer to set.
 * @param idx Location to set to.
 */
      void setPtr(Lucee::FieldPtr<T>& ptr, const int idx[NDIM]);

/**
 * Set pointer to given (i) 1D location.
 *
 * @param ptr Pointer to set.
 * @param i Location to set to.
 */
      void setPtr(Lucee::ConstFieldPtr<T>& ptr, int i) const;

/**
 * Set pointer to given (i,j) 2D location.
 *
 * @param ptr Pointer to set.
 * @param i Location to set to.
 * @param j Location to set to.
 */
      void setPtr(Lucee::ConstFieldPtr<T>& ptr, int i, int j) const;

/**
 * Set pointer to given (i,j,k) 3D location.
 *
 * @param ptr Pointer to set.
 * @param i Location to set to.
 * @param j Location to set to.
 * @param k Location to set to.
 */
      void setPtr(Lucee::ConstFieldPtr<T>& ptr, int i, int j, int k) const;

/**
 * Set pointer to given N-dimensional location.
 *
 * @param ptr Pointer to set.
 * @param idx Location to set to.
 */
      void setPtr(Lucee::ConstFieldPtr<T>& ptr, const int idx[NDIM]) const;

    private:
/**
 * Create a field attached to a given data space.
 *
 * @param rgn Region for field.
 * @param sc Start index into component.
 * @param ec End index into component.
 * @param nc Number of components in field.
 * @param lg Ghost indexes along lower index range in each dimension.
 * @param ug Ghost indexes along upper index range in each dimension.
 * @param subArr Array space to reuse.
 */
      Field(const Lucee::Region<NDIM, int>& rgn, unsigned sc, unsigned ec,
        int lg[NDIM], int ug[NDIM], Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer>& subArr);

/** Start index into components */
      unsigned scIdx;
/** Number of components */
      unsigned numComponents;
/** Region indexed by grid */
      Lucee::Region<NDIM, int> rgn;
/** Indexer over region over which field is valid */
      Lucee::RowMajorIndexer<NDIM+1> rgnIdx;
/** Lower ghost indices */
      int lowerGhost[NDIM];
/** Upper ghost indices */
      int upperGhost[NDIM];
  };

  template <unsigned NDIM, typename T>
  void
  Field<NDIM, T>::setPtr(Lucee::FieldPtr<T>& ptr, int i)
  {
    ptr.setData(&this->getRefToLoc(rgnIdx.getIndex(i,0)));
  }

  template <unsigned NDIM, typename T>
  void
  Field<NDIM, T>::setPtr(Lucee::FieldPtr<T>& ptr, int i, int j)
  {
    ptr.setData(&this->getRefToLoc(rgnIdx.getIndex(i,j,0)));
  }

  template <unsigned NDIM, typename T>
  void
  Field<NDIM, T>::setPtr(Lucee::FieldPtr<T>& ptr, int i, int j, int k)
  {
    ptr.setData(&this->getRefToLoc(rgnIdx.getIndex(i,j,k,0)));
  }

  template <unsigned NDIM, typename T>
  void
  Field<NDIM, T>::setPtr(Lucee::FieldPtr<T>& ptr, const int idx[NDIM])
  {
    int myIdx[NDIM+1];
    for (unsigned i=0; i<NDIM; ++i) myIdx[i] = idx[i];
    myIdx[NDIM] = 0;
    ptr.setData(&this->getRefToLoc(rgnIdx.getIndex(myIdx)));
  }

  template <unsigned NDIM, typename T>
  void
  Field<NDIM, T>::setPtr(Lucee::ConstFieldPtr<T>& ptr, int i) const
  {
    ptr.setData(&this->getConstRefToLoc(rgnIdx.getIndex(i,0)));
  }

  template <unsigned NDIM, typename T>
  void
  Field<NDIM, T>::setPtr(Lucee::ConstFieldPtr<T>& ptr, int i, int j) const
  {
    ptr.setData(&this->getConstRefToLoc(rgnIdx.getIndex(i,j,0)));
  }

  template <unsigned NDIM, typename T>
  void
  Field<NDIM, T>::setPtr(Lucee::ConstFieldPtr<T>& ptr, int i, int j, int k) const
  {
    ptr.setData(&this->getConstRefToLoc(rgnIdx.getIndex(i,j,k,0)));
  }

  template <unsigned NDIM, typename T>
  void
  Field<NDIM, T>::setPtr(Lucee::ConstFieldPtr<T>& ptr, const int idx[NDIM]) const
  {
    int myIdx[NDIM+1];
    for (unsigned i=0; i<NDIM; ++i) myIdx[i] = idx[i];
    myIdx[NDIM] = 0;
    ptr.setData(&this->getConstRefToLoc(rgnIdx.getIndex(myIdx)));
  }
}

#endif // LC_FIELD_H
