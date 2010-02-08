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
#include <LcColMajorIndexer.h>
#include <LcExcept.h>

namespace Lucee
{
/**
 * A field represents an array that can hold multiple values per index
 * location. Fields can be indexed directly using (i,j,k,...)
 * notation, with the last index into the field components, or using
 * field iterators.
 */
  template <unsigned NDIM, typename T>
  class Field : public Luce::Array<NDIM+1, T, Lucee::RowMajorIndexer<NDIM+1> >
  {
    public:
/**
 * Create a new field with specified shape and start indices.
 *
 * @param shape Shape of the field.
 * @param start Start indices.
 * @param nc Number of components at each index location.
 * @param initi Inital value to assigne to all components.
 */      
      Field(unsigned shape[NDIM], int start[NDIM], unsigned nc, const T& init=(T)0);
  }
}
