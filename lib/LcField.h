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
#include <LcExcept.h>
#include <LcRegion.h>
#include <LcRowMajorIndexer.h>

namespace Lucee
{
/**
 * A field represents an array that can hold multiple values per index
 * location. Fields can be indexed directly using (i,j,k,...)
 * notation, with the last index into the field components, or using
 * field iterators. Further, each field can have an extra set of
 * indices outside the main indexed region which can be used to
 * represent "ghost cells or nodes". This is useful, for example, in
 * computing finite-difference stencils on the field.
 */
  template <unsigned NDIM, typename T>
  class Field : public Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer<NDIM+1> >
  {
    public:
/**
 * Create a new field indexing given region. This constructor creates
 * an empty set of ghost indices.
 *
 * @param idxRgn Rehion indexed by array.
 * @param nc Number of components at each index location.
 * @param init Inital value to assigned to all components.
 */      
      Field(const Lucee::Region<NDIM, int>& idxRgn, unsigned nc, const T& init=(T)0);

/**
 * Create a new field indexing given region. This constructor creates
 * an empty set of ghost indices.
 *
 * @param idxRgn Rehion indexed by array.
 * @param lowerGhost Number of ghost indices along lower side of indexed region.
 * @param upperGhost Number of ghost indices along upper side of indexed region.
 * @param nc Number of components at each index location.
 * @param init Inital value to assigned to all components.
 */      
      Field(const Lucee::Region<NDIM, int>& idxRgn, 
        unsigned lowerGhost[NDIM], unsigned upperGhost[NDIM],
        unsigned nc, const T& init=(T)0);

    private:
/** Number of components */
      unsigned numComponents;
/** Lower and upper ghost indices */
      unsigned lowerGhost[NDIM], upperGhost[NDIM];
  };
}

#endif // LC_FIELD_H
