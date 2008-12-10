/**
 * @file	lckeyval.cc
 *
 * @brief	Class to hold key-value pairs.
 *
 * @version	$Id$ *
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

// lib includes
#include <lckeyval.h>

namespace Lucee
{

  KeyVal::KeyVal() 
  {
  }

  KeyVal::~KeyVal() 
  {
    // delete all entries in map
    values.erase( values.begin(), values.end() );
  }

  KeyVal::KeyVal(const KeyVal& kv)
  {
    AnyMap_t::const_iterator i;
    // copy all entries
    for (i=kv.values.begin(); i!=kv.values.end(); ++i)
      values.insert( AnyPair_t((*i).first, (*i).second) );
    // copy the type -> keys map
    copyKeys<int>(kv);
    copyKeys<double>(kv);
    copyKeys<std::string>(kv);
    copyKeys<Lucee::Any>(kv);
    copyKeys<std::vector<Lucee::Any> >(kv);
  }

  KeyVal&
  KeyVal::operator=(const KeyVal& rhs)
  {
    if (this==&rhs) return *this;

    // delete all entries in map
    values.erase( values.begin(), values.end() );
    // delete all entries in type -> key map
    Lucee::typeMapExtract<int>(typeToKeys).keys.erase(
        Lucee::typeMapExtract<int>(typeToKeys).keys.begin(),
        Lucee::typeMapExtract<int>(typeToKeys).keys.end());

    Lucee::typeMapExtract<int>(typeToKeys).keys.erase(
        Lucee::typeMapExtract<int>(typeToKeys).keys.begin(),
        Lucee::typeMapExtract<int>(typeToKeys).keys.end());

    Lucee::typeMapExtract<double>(typeToKeys).keys.erase(
        Lucee::typeMapExtract<double>(typeToKeys).keys.begin(),
        Lucee::typeMapExtract<double>(typeToKeys).keys.end());

    Lucee::typeMapExtract<std::string>(typeToKeys).keys.erase(
        Lucee::typeMapExtract<std::string>(typeToKeys).keys.begin(),
        Lucee::typeMapExtract<std::string>(typeToKeys).keys.end());

    Lucee::typeMapExtract<Lucee::Any>(typeToKeys).keys.erase(
        Lucee::typeMapExtract<Lucee::Any>(typeToKeys).keys.begin(),
        Lucee::typeMapExtract<Lucee::Any>(typeToKeys).keys.end());

    Lucee::typeMapExtract<std::vector<Lucee::Any> >(typeToKeys).keys.erase(
        Lucee::typeMapExtract<std::vector<Lucee::Any> >(typeToKeys).keys.begin(),
        Lucee::typeMapExtract<std::vector<Lucee::Any> >(typeToKeys).keys.end());

    // add entries from rhs
    AnyMap_t::const_iterator i;
    for (i=rhs.values.begin(); i!=rhs.values.end(); ++i)
      values.insert( AnyPair_t((*i).first, (*i).second) );

    // copy the type -> keys map
    copyKeys<int>(rhs);
    copyKeys<double>(rhs);
    copyKeys<std::string>(rhs);
    copyKeys<std::vector<Lucee::Any> >(rhs);

    return *this;
  }

  bool 
  KeyVal::has(const std::string& key) const
  {
    AnyMap_t::const_iterator i;
    i = values.find(key);
    return (i != values.end()) ? true : false;
  }
}
