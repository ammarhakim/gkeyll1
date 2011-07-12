/**
 * @file	LcExcept.cpp
 *
 * @brief	Class to represent exceptions in Lucee
 */

// lucee includes
#include <LcExcept.h>

namespace Lucee
{
  Except::Except() 
  {
  }

  Except::Except(const std::string& str)
    : std::exception() 
  {
    (*this) << str; 
  }

  Except::Except(const Except& ex)
  {
    (*this) << ex.exceptStrm.str();
  }

  Except::~Except() throw()
  {
  }

  Except&
  Except::operator=(const Except& ex)
  {
    if (this==&ex) return *this;
    this->exceptStrm.str(""); // zap the string
    (*this) << ex.exceptStrm.str();
    return *this;
  }

  const char*
  Except::what() const throw()
  {
    return this->exceptStrm.str().c_str();
  }
}
