/**
 * @file	LcVectorDescription.cpp
 *
 * @brief	Description of a single vector in a Lua table.
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
#include <LcVectorDescription.h>

// std includes
#include <sstream>

namespace Lucee
{

  template <typename T>
  VectorDescription<T>::VectorDescription()
    : isLengthSpecified(false), isOptnl(false), 
      varSpecified(false), var(0)
  {
  }

  template <typename T>
  VectorDescription<T>::VectorDescription(const std::vector<T>& dv)
    : isLengthSpecified(false), isOptnl(true), defVector(dv), 
      varSpecified(false), var(0)
  {
  }

  template <typename T>
  void
  VectorDescription<T>::fillVarWithVector(const std::vector<T>& vec)
  {
    var->clear(); // clear in case called twice
    for (unsigned i=0; i<vec.size(); ++i)
      var->push_back(vec[i]);
  }

  template <typename T>
  void
  VectorDescription<T>::fillVarWithOptional()
  {
    var->clear(); // clear in case called twice
    for (unsigned i=0; i<defVector.size(); ++i)
      var->push_back(defVector[i]);
  }

  template <typename T>
  std::pair<bool, std::string>
  VectorDescription<T>::checkVector(const std::vector<T>& vec)
  {
    bool pass = true;
    std::ostringstream errMsg;
// check for size of vector
    if (isLengthSpecified)
    { // check length
      if (vec.size() != length)
      {
        pass = false;
        errMsg << "** (Vector not of correct size. Should have length "
               << length << ". Supplied vector is of length " << vec.size() << ")";
      }
    }
// check each element in vector
    for (unsigned i=0; i<vec.size(); ++i)
    {
      std::pair<bool, std::string> ss = lastValDescr.checkValue(vec[i]);
      if (ss.first == false)
      {
        pass = false;
        errMsg << "** (For element " << i << " " << ss.second << ")" << std::endl;
      }
    }

    return std::pair<bool, std::string>(pass, errMsg.str());
  }

  template <typename T>
  VectorDescription<T>&
  VectorDescription<T>::setLength(unsigned sz)
  {
    isLengthSpecified = true;
    length = sz;
    return *this;
  }

  template <typename T>
  VectorDescription<T>&
  VectorDescription<T>::setHelp(const std::string& hlp)
  {
    lastValDescr.setHelp(hlp);
    return *this;
  }

  template <typename T>
  VectorDescription<T>&
  VectorDescription<T>::setMinValue(const T& mv)
  {
    lastValDescr.setMinValue(mv);
    return *this;
  }

  template <typename T>
  VectorDescription<T>&
  VectorDescription<T>::setMaxValue(const T& mv)
  {
    lastValDescr.setMaxValue(mv);
    return *this;
  }

  template <typename T>
  VectorDescription<T>&
  VectorDescription<T>::setOneOf(const std::vector<T>& onef)
  {
  }

  template <typename T>
  VectorDescription<T>&
  VectorDescription<T>::setVar(std::vector<T>* v)
  {
    varSpecified = true;
    var = v;
  }

// instantiations
  template class VectorDescription<int>;
  template class VectorDescription<double>;
  template class VectorDescription<std::string>;
}
