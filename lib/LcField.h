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
  class Field : public Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer<NDIM+1> >
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
 * Create a new pointer object to elements in field.
 */
      Lucee::FieldPtr<T> createPtr();

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

    private:
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
  Field<NDIM, T>::Field(const Lucee::Region<NDIM, int>& rgn, unsigned nc, const T& init)
    : Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer<NDIM+1> >(rgn.inflate(0, nc), init),
      numComponents(nc), rgn(rgn), rgnIdx(rgn.inflate(0, nc))
  {
    for (unsigned i=0; i<NDIM; ++i)
    {
      lowerGhost[i] = 0;
      upperGhost[i] = 0;
    }
  }

  template <unsigned NDIM, typename T>
  Field<NDIM, T>&
  Field<NDIM, T>::operator=(const T& val)
  {
// simply call base class assignment operator    
    Array<NDIM+1, T, Lucee::RowMajorIndexer<NDIM+1> >::operator=(val);
    return *this;
  }

  template <unsigned NDIM, typename T>
  Lucee::FieldPtr<T>
  Field<NDIM, T>::createPtr()
  {
    int start[NDIM];
    for (unsigned i=0; i<NDIM; ++i)
      start[i] = rgn.getLower(i);
    unsigned loc = rgnIdx.getGenLowIndex(start);
    return Lucee::FieldPtr<T>(numComponents, &this->getRefToLoc(loc));
  }

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
}

#endif // LC_FIELD_H
