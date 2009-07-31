/**
 * @file	LcObject.cc
 *
 * @brief	Base class for all Lucee objects
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// lucee includes
#include <LcObject.h>

namespace Lucee
{
  Object::Object(const std::string& name)
    : name(name) 
  {
  }

  Object::~Object()
  {
  }

  std::string Object::getName() const 
  {
    return name;
  }
}
