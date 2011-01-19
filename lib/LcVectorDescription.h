/**
 * @file	LcVectorDescription.h
 *
 * @brief	Description of a single vector in a Lua table.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_VECTOR_DESCRIPTION_H
#define LC_VECTOR_DESCRIPTION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcValueDescription.h>

// std includes
#include <string>
#include <vector>

namespace Lucee
{
/**
  Description of a single vector in a Lua table. All elements in the
  vector must have the same type. This description is part of the Lua
  script validation file.
 */
  template <typename T>
  class VectorDescription
  {
    public:
/**
 * Create a new vector description object.
 */
      VectorDescription();

    private:
/** Name of vector */
      std::string name;
/** Description of each element in vector */
      std::vector<Lucee::ValueDescription<T> > valDescr;
/** Description of final element (in case of variable size vector) */
      Lucee::ValueDescription<T> lastValDescr;
  };
}

#endif // LC_VECTOR_DESCRIPTION_H
