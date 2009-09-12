/**
 * @file	LcKeyVal.cc
 *
 * @brief	Class to hold key-value pairs.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee include
#include <LcKeyVal.h>

namespace Lucee
{
  KeyVal::KeyVal()
  {
  }

  KeyVal::~KeyVal()
  {
  }

  KeyVal::KeyVal(const KeyVal& kv)
  {
  }

  KeyVal&
  KeyVal::operator=(const KeyVal& kv)
  {
    return *this;
  }
}
