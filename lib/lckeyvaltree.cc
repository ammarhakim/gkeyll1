/**
 * @file	lckeyvaltree.cc
 *
 * @brief	Class to hold a tree of key-value pairs.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

// lib includes
#include <lckeyvaltree.h>
#include <lcexcept.h>

namespace Lucee
{

  KeyValTree::KeyValTree()
    : Lucee::KeyVal(), name("")
  {}

  KeyValTree::KeyValTree(const std::string& name)
    : Lucee::KeyVal(), name(name)
  {}

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

    for (i=kvt.kvTreeMap.begin(); i!=kvt.kvTreeMap.end(); ++i)
    {
      KeyValTree *kvtp = new KeyValTree(*(i->second));
      kvTreeMap.insert( KeyValTreePair_t( i->first, kvtp ) );
      // check if cryptset has Type attribute
      if (kvtp->has("Type"))
      {
        std::string type = kvtp->get<std::string>("Type");
        // add name to type map
        TypeMap_t::iterator i;
        i = typeMap.find(type);
        if (i != typeMap.end())
          i->second.push_back(kvtp->getName());
        else
        {
          std::vector<std::string> name;
          name.push_back(kvtp->getName());
          typeMap.insert( TypePair_t(type, name) );
        }
      }
    }
  }

  KeyValTree&
  KeyValTree::operator=(const KeyValTree& kvt)
  {
    if (this==&kvt) return *this;

    // call base class assignment operator
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
      // check if cryptset has Type attribute
      if (kvtp->has("Type"))
      {
        std::string type = kvtp->get<std::string>("Type");
        // add name to type map
        TypeMap_t::iterator i;
        i = typeMap.find(type);
        if (i != typeMap.end())
          (i->second).push_back(kvtp->getName());
        else
        {
          std::vector<std::string> name;
          name.push_back(kvtp->getName());
          typeMap.insert( TypePair_t(type, name) );
        }
      }
    }

    return *this;
  }
  
  std::string
  KeyValTree::getName() const
  {
    return name;
  }

  void 
  KeyValTree::addSet(const KeyValTree& kvt) 
  {
    // make a deep copy of the set
    KeyValTree *kvtp = new KeyValTree(kvt);
    kvTreeMap.insert( KeyValTreePair_t(kvt.getName(), kvtp) );
    // check if cryptset has Type attribute
    if (kvt.has("Type"))
    {
      std::string type = kvt.get<std::string>("Type");
      // add name to type map
      TypeMap_t::iterator i;
      i = typeMap.find(type);
      if (i != typeMap.end())
        (i->second).push_back(kvt.getName());
      else
      {
        std::vector<std::string> name;
        name.push_back(kvt.getName());
        typeMap.insert( TypePair_t(type, name) );
      }
    }
  }

  bool
  KeyValTree::hasSet(const std::string& name) const 
  {
    // retrieve the set with given name
    KeyValTreeMap_t::const_iterator i;
    i = kvTreeMap.find(name);
    return (i != kvTreeMap.end()) ? true : false;
  }

  const 
  KeyValTree& 
  KeyValTree::getSet(const std::string& name) const 
  {
    // retrieve the set with given name
    KeyValTreeMap_t::const_iterator i;
    i = kvTreeMap.find(name);
    if (i != kvTreeMap.end())
      return *(i->second);
    Lucee::Except ex;
    ex << "KeyValTree set " << name << " not found";
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
