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

    private:
/** Number of components */
      unsigned numComponents;
/** Region indexed by grid */
      Lucee::Region<NDIM, int> rgn;
/** Indexer over region over which field is valid */
      Lucee::RowMajorIndexer<NDIM> rgnIdx;
/** Lower and upper ghost indices */
      unsigned lowerGhost[NDIM], upperGhost[NDIM];
  };
}

#endif // LC_FIELD_H
