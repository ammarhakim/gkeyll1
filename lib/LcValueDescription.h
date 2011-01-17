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
 * Create a new value with given name.
 *
 * @param nm Name of value as it should appear in Lua script.
 */
      ValueDescription(const std::string& nm);

/**
 * Create a new value with given name. This constructor create a value
 * that is option and hence needs a default value.
 *
 * @param nm Name of value as it should appear in Lua script.
 */
      ValueDescription(const std::string& nm, const T& dv);

/**
 * Set help string.
 *
 * @param hlp Help string.
 */
      void setHelp(const std::string& hlp);

/**
 * Set minimum possible value (inclusive).
 *
 * @param mv Minimum value.
 */
      void setMinValue(const T& mv);

/**
 * Set maximum possible value (inclusive).
 *
 * @param mv Maximum value.
 */
      void setMaxValue(const T& mv);

/**
 * Set possible values that this can take. This takes precedence over
 * min/max values, i.e. if list of possible values is specified, then
 * min/max checkes are skipped.
 *
 * @param onef Value must be one of these.
 */
      void setOneOf(const std::vector<T>& onef);

    private:
/** Name of value as it appears in the Lua script */
      std::string name;
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
