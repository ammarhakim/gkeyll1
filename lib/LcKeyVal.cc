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
    : mapCntr(new DataStoreMapCntr)
  {
  }

  KeyVal::KeyVal(const KeyVal& kv)
  { // shallow copy
    this->mapCntr = kv.mapCntr;
  }

  KeyVal&
  KeyVal::operator=(const KeyVal& kv)
  { // shallow copy
    if (this == &kv)
      return *this;
    this->mapCntr = kv.mapCntr;
    return *this;
  }

  KeyVal
  KeyVal::duplicate() const
  {
    KeyVal kv;
// copy stuff from us into kv
    copyFromSet<int>(kv, *this);
    copyFromSet<double>(kv, *this);
    copyFromSet<std::string>(kv, *this);

    copyFromSet<std::vector<int> >(kv, *this);
    copyFromSet<std::vector<double> >(kv, *this);
    copyFromSet<std::vector<std::string> >(kv, *this);

    return kv;
  }
}
