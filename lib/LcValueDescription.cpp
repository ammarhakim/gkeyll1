/**
 * @file	LcValueDescription.cpp
 *
 * @brief	Description of a single value in a Lua table.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcValueDescription.h>

namespace Lucee
{
  template <typename T>
  ValueDescription<T>::ValueDescription(const std::string& nm)
    : name(nm), isOptional(false), varSpecified(false), var(0),
      isOneOf(false), isMinSet(false), isMaxSet(false)
  {
  }

  template <typename T>
  ValueDescription<T>::ValueDescription(const std::string& nm, const T& dv)
    : name(nm), isOptional(true), defValue(dv), varSpecified(false), var(0),
      isOneOf(false), isMinSet(false), isMaxSet(false)
  {
  }

  template <typename T>
  void
  ValueDescription<T>::setHelp(const std::string& hlp)
  {
    help = hlp;
  }

  template <typename T>
  void
  ValueDescription<T>::setMinValue(const T& mv)
  {
    isMinSet = true;
    minVal = mv;
  }
  
  template <typename T>
  void
  ValueDescription<T>::setMaxValue(const T& mv)
  {
    isMaxSet = true;
    maxVal = mv;
  }

  template <typename T>
  void
  ValueDescription<T>::setOneOf(const std::vector<T>& onef)
  {
    isOneOf = true;
    oneOf.clear(); // clear in case this is called multiple times
    oneOf = onef;
  }

// instantiations
  template class ValueDescription<int>;
  template class ValueDescription<float>;
  template class ValueDescription<std::string>;
}
