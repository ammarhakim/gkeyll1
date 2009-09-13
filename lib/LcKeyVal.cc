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

  KeyVal::KeyVal(const KeyVal& kv)
  {
    copyFromSet<int>(kv);
    copyFromSet<double>(kv);
    copyFromSet<std::string>(kv);
    copyFromSet<std::vector<int> >(kv);
    copyFromSet<std::vector<double> >(kv);
    copyFromSet<std::vector<std::string> >(kv);
  }

  KeyVal&
  KeyVal::operator=(const KeyVal& kv)
  {
    if (this == &kv)
      return *this;
// clear each of our sets
    Loki::Field<int>(dataStoreTypeMap).dataMap.clear();
    Loki::Field<double>(dataStoreTypeMap).dataMap.clear();
    Loki::Field<std::string>(dataStoreTypeMap).dataMap.clear();
    Loki::Field<std::vector<int> >(dataStoreTypeMap).dataMap.clear();
    Loki::Field<std::vector<double> >(dataStoreTypeMap).dataMap.clear();
    Loki::Field<std::vector<std::string> >(dataStoreTypeMap).dataMap.clear();
// now copy stuff over
    copyFromSet<int>(kv);
    copyFromSet<double>(kv);
    copyFromSet<std::string>(kv);
    copyFromSet<std::vector<int> >(kv);
    copyFromSet<std::vector<double> >(kv);
    copyFromSet<std::vector<std::string> >(kv);

    return *this;
  }
}
