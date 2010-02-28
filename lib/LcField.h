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
 * Region indexed by field.
 *
 * @return region indexed by field.
 */
      Lucee::Region<NDIM, int> getRegion() const { return rgn; }

/**
 * Get a lower-dimensional field object. The returned object shares
 * data with this object. The return field has access to all
 * components stored in this field.
 *
 * @param defDims Dimensions to remove.
 * @param defDimsIdx Index of removed dimensions.
 * @return Lower dimensional field object.
 */
      template <unsigned RDIM, T>
      Field<RDIM, T>
      getLowDimField(const unsigned defDims[NDIM-RDIM], const int defDimsIdx[NDIM-RDIM])
      {
      }

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
 * Get a field with same indexed region as this one, but with 
 */
      Field<NDIM, T> getSubCompField(unsigned sc, unsigned ec);

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
 * @param nc Number of components in field.
 * @param subArr Array space to reuse.
 */
      Field(const Lucee::Region<NDIM, int>& rgn, unsigned nc,
        Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer>& subArr);

/** Number of components */
      unsigned numComponents;
/** Region indexed by grid */
      Lucee::Region<NDIM, int> rgn;
/** Indexer over region over which field is valid */
      Lucee::RowMajorIndexer<NDIM+1> rgnIdx;
/** Lower and upper ghost indices */
      unsigned lowerGhost[NDIM], upperGhost[NDIM];
  };

  template <unsigned NDIM, typename T>
  void
  Field<NDIM, T>::setPtr(Lucee::FieldPtr<T>& ptr, int i)
  {
    ptr.setData(&this->getRefToLoc(rgnIdx.getLowIndex(i)));
  }

  template <unsigned NDIM, typename T>
  void
  Field<NDIM, T>::setPtr(Lucee::FieldPtr<T>& ptr, int i, int j)
  {
    ptr.setData(&this->getRefToLoc(rgnIdx.getLowIndex(i,j)));
  }

  template <unsigned NDIM, typename T>
  void
  Field<NDIM, T>::setPtr(Lucee::FieldPtr<T>& ptr, int i, int j, int k)
  {
    ptr.setData(&this->getRefToLoc(rgnIdx.getLowIndex(i,j,k)));
  }

  template <unsigned NDIM, typename T>
  void
  Field<NDIM, T>::setPtr(Lucee::FieldPtr<T>& ptr, const int idx[NDIM])
  {
    ptr.setData(&this->getRefToLoc(rgnIdx.getGenLowIndex(idx)));
  }

  template <unsigned NDIM, typename T>
  void
  Field<NDIM, T>::setPtr(Lucee::ConstFieldPtr<T>& ptr, int i) const
  {
    ptr.setData(&this->getConstRefToLoc(rgnIdx.getLowIndex(i)));
  }

  template <unsigned NDIM, typename T>
  void
  Field<NDIM, T>::setPtr(Lucee::ConstFieldPtr<T>& ptr, int i, int j) const
  {
    ptr.setData(&this->getConstRefToLoc(rgnIdx.getLowIndex(i,j)));
  }

  template <unsigned NDIM, typename T>
  void
  Field<NDIM, T>::setPtr(Lucee::ConstFieldPtr<T>& ptr, int i, int j, int k) const
  {
    ptr.setData(&this->getConstRefToLoc(rgnIdx.getLowIndex(i,j,k)));
  }

  template <unsigned NDIM, typename T>
  void
  Field<NDIM, T>::setPtr(Lucee::ConstFieldPtr<T>& ptr, const int idx[NDIM]) const
  {
    ptr.setData(&this->getConstRefToLoc(rgnIdx.getGenLowIndex(idx)));
  }
}

#endif // LC_FIELD_H
