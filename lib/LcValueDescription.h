/**
 * @file	LcValueDescription.h
 *
 * @brief	Description of a single value in a Lua table.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_VALUE_DESCRIPTION_H
#define LC_VALUE_DESCRIPTION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// std includes
#include <string>
#include <vector>

namespace Lucee
{
/**
 * Description of a single value in a Lua table. This description is
 * part of the Lua script validation file.
 */
  template <typename T>
  class ValueDescription
  {
    public:
/**
 * Create a new value description object.
 */
      ValueDescription();

/**
 * Create a new value description object. This constructor creates a
 * value that is option and hence needs a default value.
 *
 * @param dv Default value.
 */
      ValueDescription(const T& dv);

/**
 * Set help string.
 *
 * @param hlp Help string.
 * @return reference to this object.
 */
      ValueDescription<T>& setHelp(const std::string& hlp);

/**
 * Set minimum possible value (inclusive).
 *
 * @param mv Minimum value.
 * @return reference to this object.
 */
      ValueDescription<T>& setMinValue(const T& mv);

/**
 * Set maximum possible value (inclusive).
 *
 * @param mv Maximum value.
 * @return reference to this object.
 */
      ValueDescription<T>& setMaxValue(const T& mv);

/**
 * Set possible values that this can take. This takes precedence over
 * min/max values, i.e. if list of possible values is specified, then
 * min/max checkes are skipped.
 *
 * @param onef Value must be one of these.
 * @return reference to this object.
 */
      ValueDescription<T>& setOneOf(const std::vector<T>& onef);

/**
 * Set pointer to variable that will be set.
 *
 * @param var pointer to dat that should be set.
 * @return reference to this object.
 */
      ValueDescription<T>& setVar(T* var);

    private:
/** Help string for value */
      std::string help;
/** Is this value optional? */
      bool isOptional;
/** Default value if it is optional */
      T defValue;
/** Was a settable variable specified? */
      bool varSpecified;
/** Pointer to settable variable */
      T *var;
/** Was oneOf set? */
      bool isOneOf;
/** Can be one of these values */
      std::vector<T> oneOf;
/** Was min set? */
      bool isMinSet;
/** Minimum value */
      T minVal;
/** Was max set? */
      bool isMaxSet;
/** Minimum value */
      T maxVal;
  };
}

#endif // LC_VALUE_DESCRIPTION_H
