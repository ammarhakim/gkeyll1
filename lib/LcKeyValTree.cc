/**
 * @file	LcKeyValTree.cc
 *
 * @brief	Class to hold a tree of key-value pairs.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// lib includes
#include <LcKeyValTree.h>
#include <LcExcept.h>

namespace Lucee
{
  KeyValTree::KeyValTree(const std::string& name, const std::string& type,
    const std::string& kind)
    : Lucee::KeyVal(), name(name), type(type), kind(kind)
  {
  }

  KeyValTree::~KeyValTree()
  {
    KeyValTreeMap_t::iterator i;
    for (i=kvTreeMap.begin(); i!=kvTreeMap.end(); ++i)
      delete i->second;
    kvTreeMap.erase( kvTreeMap.begin(), kvTreeMap.end() );
  }

  KeyValTree::KeyValTree(const KeyValTree& kvt)
    : Lucee::KeyVal(kvt)
  {
    KeyValTreeMap_t::const_iterator i;
    name = kvt.name;
    type = kvt.type;
    kind = kvt.kind;

    for (i=kvt.kvTreeMap.begin(); i!=kvt.kvTreeMap.end(); ++i)
    {
      KeyValTree *kvtp = new KeyValTree(*(i->second));
      kvTreeMap.insert( KeyValTreePair_t( i->first, kvtp ) );
    }
  }

  KeyValTree&
  KeyValTree::operator=(const KeyValTree& kvt)
  {
    if (this==&kvt) return *this;

    Lucee::KeyVal::operator= (kvt);

    KeyValTreeMap_t::const_iterator i;
// delete all entries in map
    for (i=kvTreeMap.begin(); i!=kvTreeMap.end(); ++i)
      delete i->second;
    kvTreeMap.erase( kvTreeMap.begin(), kvTreeMap.end() );

// copy entries from kvt
    name = kvt.name;
    for (i=kvt.kvTreeMap.begin(); i!=kvt.kvTreeMap.end(); ++i)
    {
      KeyValTree *kvtp = new KeyValTree(*(i->second));
      kvTreeMap.insert( KeyValTreePair_t( i->first, kvtp ) );
    }

    return *this;
  }
  
  bool
  KeyValTree::addTree(const KeyValTree& kvt) 
  {
    KeyValTree *kvtp = new KeyValTree(kvt);
    return kvTreeMap.insert(
      KeyValTreePair_t(kvt.getName(), kvtp)).second;
  }

  bool
  KeyValTree::hasTree(const std::string& name) const 
  {
    KeyValTreeMap_t::const_iterator i = kvTreeMap.find(name);
    return (i != kvTreeMap.end()) ? true : false;
  }

  const 
  KeyValTree& 
  KeyValTree::getTree(const std::string& name) const 
  {
    KeyValTreeMap_t::const_iterator i;
    i = kvTreeMap.find(name);
    if (i != kvTreeMap.end())
      return *(i->second);
    Lucee::Except ex;
    ex << "KeyValTree::getTree: set " << name << " not found";
    throw ex;
  }

  std::vector<std::string> 
  KeyValTree::getNamesOfType(const std::string& type) const
  {
    TypeMap_t::const_iterator i;
    i = typeMap.find(type);
    if (i != typeMap.end())
      return i->second;
    else
      return std::vector<std::string>();
  }
}
