/**
 * @file	lcexcept.cc
 *
 * @brief	Class to represent exceptions in Lucee
 *
 * @version	$Id$ *
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

// lib includes
#include <lcexcept.h>

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
    : std::exception(ex) 
  {
    (*this) << ex.str();
  }

  Except::~Except() throw()
  {
  }

  Except& Except::operator=(const Except& ex)
  {
    if (this==&ex) return *this;
    (*this).str(""); // zap the string
    (*this) << ex.str();
    return *this;
  }
}
