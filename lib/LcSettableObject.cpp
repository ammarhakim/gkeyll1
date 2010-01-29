/**
 * @file	LcSettableObject.cc
 *
 * @brief	Base class for all Lucee objects which can be set from name/values pairs.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// lucee includes
#include <LcExcept.h>
#include <LcSettableObject.h>

namespace Lucee
{
  SettableObject::SettableObject(const std::string& name)
    : name(name) 
  {
  }

  SettableObject::~SettableObject()
  {
  }

  std::string
  SettableObject::getName() const 
  {
    return name;
  }
}
