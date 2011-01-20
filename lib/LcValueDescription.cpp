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

// std includes
#include <sstream>

namespace Lucee
{
  template <typename T>
  ValueDescription<T>::ValueDescription()
    : isOptnl(false), varSpecified(false), var(0),
      isOneOf(false), isMinSet(false), isMaxSet(false)
  {
  }

  template <typename T>
  ValueDescription<T>::ValueDescription(const T& dv)
    : isOptnl(true), defValue(dv), varSpecified(false), var(0),
      isOneOf(false), isMinSet(false), isMaxSet(false)
  {
  }

  template <typename T>
  void
  ValueDescription<T>::fillVarWithValue(const T& val)
  {
// set if specified
    if (varSpecified) *var = val;
  }

  template <typename T>
  void
  ValueDescription<T>::fillVarWithOptional()
  {
    if (isOptnl)
// set if specified
      if (varSpecified) *var = defValue;
  }

  template <typename T>
  std::pair<bool, std::string>
  ValueDescription<T>::checkValue(const T& val)
  {
    bool pass = true;
    std::ostringstream errMsg;
    if (isOneOf)
    { // check if value if one of specified list
      bool isOneOfPass = false;
      for (int i=0; i<oneOf.size(); ++i)
        if (val == oneOf[i])
          isOneOfPass = true;
      if (isOneOfPass == false)
      {
        pass = false;
        errMsg << "Value '" << val << "' is not one of [";
        for (int i=0; i<oneOf.size()-1; ++i)
          errMsg << oneOf[i] << ", ";
        errMsg << oneOf[oneOf.size()-1] << "]";
      }
    }
    else
    {
      if (isMinSet)
      {
        if (val < minVal)
        {
          pass = false;
          errMsg << "Value '" << val << "' is less than " << minVal;
        }
      }
      if (isMaxSet)
      {
        if (val > maxVal)
        {
          pass = false;
          errMsg << "Value '" << val << "' is greater than " << maxVal;
        }
      }
    }

    return std::pair<bool, std::string>(pass, errMsg.str());
  }

  template <typename T>
  ValueDescription<T>&
  ValueDescription<T>::setHelp(const std::string& hlp)
  {
    help = hlp;
    return *this;
  }

  template <typename T>
  ValueDescription<T>&
  ValueDescription<T>::setMinValue(const T& mv)
  {
    isMinSet = true;
    minVal = mv;
    return *this;
  }
  
  template <typename T>
  ValueDescription<T>&
  ValueDescription<T>::setMaxValue(const T& mv)
  {
    isMaxSet = true;
    maxVal = mv;
    return *this;
  }

  template <typename T>
  ValueDescription<T>&
  ValueDescription<T>::setOneOf(const std::vector<T>& onef)
  {
    isOneOf = true;
    oneOf.clear(); // clear in case this is called multiple times
    oneOf = onef;
    return *this;
  }

  template <typename T>
  ValueDescription<T>&
  ValueDescription<T>::setVar(T* v)
  {
    varSpecified = true;
    var = v;
  }

// instantiations
  template class ValueDescription<int>;
  template class ValueDescription<double>;
  template class ValueDescription<std::string>;
}
