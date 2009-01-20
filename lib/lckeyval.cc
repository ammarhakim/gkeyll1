/**
 * @file	lckeyval.cc
 *
 * @brief	Class to hold key-value pairs.
 *
 * @version	$Id$
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
    copyKeys<float>(kv);
    copyKeys<double>(kv);
    copyKeys<std::string>(kv);

    copyKeys<std::vector<int> >(kv);
    copyKeys<std::vector<float> >(kv);
    copyKeys<std::vector<double> >(kv);
    copyKeys<std::vector<std::string> >(kv);
  }

  KeyVal&
  KeyVal::operator=(const KeyVal& kv)
  {
    if (this==&kv) return *this;

    // delete all entries in map
    values.erase( values.begin(), values.end() );
    // delete all entries in type -> key map
    Loki::Field<int>(typeToKeys).keys.erase(
        Loki::Field<int>(typeToKeys).keys.begin(),
        Loki::Field<int>(typeToKeys).keys.end());

    Loki::Field<int>(typeToKeys).keys.erase(
        Loki::Field<int>(typeToKeys).keys.begin(),
        Loki::Field<int>(typeToKeys).keys.end());

    Loki::Field<double>(typeToKeys).keys.erase(
        Loki::Field<double>(typeToKeys).keys.begin(),
        Loki::Field<double>(typeToKeys).keys.end());

    Loki::Field<std::string>(typeToKeys).keys.erase(
        Loki::Field<std::string>(typeToKeys).keys.begin(),
        Loki::Field<std::string>(typeToKeys).keys.end());

    Loki::Field<std::vector<int> >(typeToKeys).keys.erase(
        Loki::Field<std::vector<int> >(typeToKeys).keys.begin(),
        Loki::Field<std::vector<int> >(typeToKeys).keys.end());

    Loki::Field<std::vector<float> >(typeToKeys).keys.erase(
        Loki::Field<std::vector<float> >(typeToKeys).keys.begin(),
        Loki::Field<std::vector<float> >(typeToKeys).keys.end());

    Loki::Field<std::vector<double> >(typeToKeys).keys.erase(
        Loki::Field<std::vector<double> >(typeToKeys).keys.begin(),
        Loki::Field<std::vector<double> >(typeToKeys).keys.end());

    Loki::Field<std::vector<std::string> >(typeToKeys).keys.erase(
        Loki::Field<std::vector<std::string> >(typeToKeys).keys.begin(),
        Loki::Field<std::vector<std::string> >(typeToKeys).keys.end());

    // add entries from kv
    AnyMap_t::const_iterator i;
    for (i=kv.values.begin(); i!=kv.values.end(); ++i)
      values.insert( AnyPair_t((*i).first, (*i).second) );

    // copy the type -> keys map
    copyKeys<int>(kv);
    copyKeys<float>(kv);
    copyKeys<double>(kv);
    copyKeys<std::string>(kv);

    copyKeys<std::vector<int> >(kv);
    copyKeys<std::vector<float> >(kv);
    copyKeys<std::vector<double> >(kv);
    copyKeys<std::vector<std::string> >(kv);

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
