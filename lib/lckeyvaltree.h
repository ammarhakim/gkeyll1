/**
 * @file	lckeyvaltree.h
 *
 * @brief	Class to hold a tree of key-value pairs.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_KEYVAL_TREE_H
#define LC_KEYVAL_TREE_H

// lib includes
#include <lckeyval.h>

// std incudes
#include <map>
#include <string>
#include <vector>

namespace Lucee
{
/**
 * Lucee::KeyValTree extends Lucee::KeyVal by providing, in addition
 * to key-value pairs, an set of named Lucee::KeyValTrees, thus
 * providing a powerful way of representing hierarchical data. In
 * effect the Lucee::KeyValTree is a n-ary tree with non-unique keys.
 *
 * Copying (or assignment of) a Lucee::KeyValTree can be expensive as the
 * complete tree is copied (deep-copy semantics). In general this is
 * not a problem as the creation of the set is usually done at
 * start-up time and the resulting set is never modified.
 *
 * Like its parent Lucee::KeyVal, Lucee::KeyValTree is also immutable:
 * sets once added can not be deleted or modified.
 */
  class KeyValTree : public Lucee::KeyVal
  {
    public:
/**
 * Create empty key-value tree with given name.
 *
 * @param name Name of the set
 */
      KeyValTree(const std::string& name="KeyValTree");

/**
 * Destroy key-value tree object
 */
      virtual ~KeyValTree();

/**
 * Makes a deep copy of the supplied tree
 *
 * @param kvt Key-value tree to copy
 */
      KeyValTree(const KeyValTree& kvt);

/**
 * Assignment operator: makes a deep copy of supplied tree
 *
 * @param kvt Key-value tree to copy
 * @return reference to copy
 */
      KeyValTree& operator=(const KeyValTree& kvt);

/**
 * Set name of tree
 *
 * @param name Name of tree
 */
      void setName(const std::string& name);

/**
 * Name of tree
 *
 * @return Name of tree
 */
      std::string getName() const;

/**
 * Add a new key-value tree to this tree.
 *
 * @param kvt Add given key-value tree to this set.
 */
      void addSet(const KeyValTree& kvt);

/**
 * Check if tree with given name exists in this tree.
 *
 * @param name Name of tree to test for.
 * @return true if tree exists, false otherwise.
 */
      bool hasSet(const std::string& name) const;

/**
 * Get key-value tree with given name. The returned set is immutable.
 *
 * @param name Name of key-value tree to get.
 * @return immutable reference to tree.
 */
      const KeyValTree& getSet(const std::string& name) const;

/**
 * Return list of key-value tree names with the given type. The "type"
 * of a tree is define by the value of the key "Type" in the key-value
 * pairs. Not all trees will have a type.
 *
 * @param type Type name
 * @return names of all trees with given type
 */
      std::vector<std::string> getNamesOfType(const std::string& type) const;

/** Typedef for string->Lucee::KeyValTree map */
      typedef std::map<std::string, Lucee::KeyValTree*, std::less<std::string> > KeyValTreeMap_t;
/** Typedef for string->Lucee::KeyValTree pair */
      typedef std::pair<std::string, Lucee::KeyValTree*> KeyValTreePair_t;

/** Typedef for type -> list of Lucee::KeyValTrees map */
      typedef std::map<std::string, std::vector<std::string> > TypeMap_t;
/** Typedef of type -> list of Lucee::KeyValTrees pair */
      typedef std::pair<std::string, std::vector<std::string> > TypePair_t;

    private:
/** Name of tree */
      std::string name;
/** Map of tree names to trees */
      KeyValTreeMap_t kvTreeMap;
/** Map of tree-type to tree names */
      TypeMap_t typeMap;
  };
}

#endif // LC_KEYVAL_TREE_H
