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

/**
 * Create a new optional vector description object.
 *
 * @param dv Default value of vector.
 */
      VectorDescription(const std::vector<T>& dv);

/**
 * Set expected size of vector.
 *
 * @param sz Size of vector.
 * @return reference to this object.
 */
      VectorDescription<T>& setLength(unsigned sz);

/**
 * Set help string.
 *
 * @param hlp Help string.
 * @return reference to this object.
 */
      VectorDescription<T>& setHelp(const std::string& hlp);

/**
 * Set minimum possible value (inclusive).
 *
 * @param mv Minimum value.
 * @return reference to this object.
 */
      VectorDescription<T>& setMinValue(const T& mv);

/**
 * Set maximum possible value (inclusive).
 *
 * @param mv Maximum value.
 * @return reference to this object.
 */
      VectorDescription<T>& setMaxValue(const T& mv);

/**
 * Set possible values that this can take. This takes precedence over
 * min/max values, i.e. if list of possible values is specified, then
 * min/max checkes are skipped.
 *
 * @param onef Value must be one of these.
 * @return reference to this object.
 */
      VectorDescription<T>& setOneOf(const std::vector<T>& onef);

/**
 * Set pointer to variable that will be set.
 *
 * @param var pointer to dat that should be set.
 * @return reference to this object.
 */
      VectorDescription<T>& setVar(std::vector<T>* var);

    private:
/** Is length of vector specified? */
      bool isLengthSpecified;
/** Length of vector */
      unsigned length;
/** Is this vector optional? */
      bool isOptional;
/** Default value if it is optional */
      std::vector<T> defValue;
/** Was a settable variable specified? */
      bool varSpecified;
/** Pointer to settable variable */
      std::vector<T> *var;
/** Description of all elements in vector (for now I am assuming all
 * elements have same discription) */
      Lucee::ValueDescription<T> lastValDescr;
  };
}

#endif // LC_VECTOR_DESCRIPTION_H
