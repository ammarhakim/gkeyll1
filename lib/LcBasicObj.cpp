/**
 * @file	LcBasicObj.cpp
 *
 * @brief	Interface class for basic Lucee objects.
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
#include <LcBasicObj.h>

// std includes
#include <string>

namespace Lucee
{
  BasicObj::BasicObj(const std::string& nm)
    : nm(nm) 
  {
  }

  BasicObj::~BasicObj() 
  {
  }

  std::string
  BasicObj::getName() const { return nm; }

  void
  BasicObj::setName(const std::string& name) 
  {
    nm = name;
  }
}
