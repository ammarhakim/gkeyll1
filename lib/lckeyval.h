/**
 * @file	lckeyval.h
 *
 * @brief	Class to hold key-value pairs.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_KEYVAL_H
#define LC_KEYVAL_H

// lib includes
#include <lcany.h>
#include <lcdatatypes.h>
#include <lcexcept.h>
#include <lctypelist.h>

// std includes
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace Lucee
{
/**
 * Lucee::KeyVal provides a container to store/retrive key-value
 * pairs. The key are strings, while the values can be of arbitrary
 * types supporting the copy constructor. For reasons of efficiency
 * Lucee::KeyVal should be used to store only "light-weight" objects
 * like native-type objects (int, double, string etc.) or small
 * objects which are not expensive to copy.
 *
 * Lucee::KeyVal objects are immutable, i.e, values once added can not
 * be removed or modifed.
 */
  class KeyVal
  {
    public:
/**
 * Create a new KeyVal object.
 */      
      KeyVal();

/**
 * Destroy the object
 */
      virtual ~KeyVal();

/**
 * Copy ctor: make a deep copy of supplied KeyVal class
 *
 * @param kv KeyVal object to copy.
 */
      KeyVal(const KeyVal& kv);

/**
 * Assignment operator: makes a deep copy.
 *
 * @param kv KeyVal object to copy.
 * @return reference to assigned object.
 */
      KeyVal& operator=(const KeyVal& kv);

/**
 * Insert a key, value pair. Returns true if the insertion worked,
 * false otherwise
 *
 * @param key Key of object.
 * @param value Value of object.
 * @return true, if insertion worked, false otherwise.
 */
      template<typename VALUETYPE>
      bool add(const std::string& key, VALUETYPE value) {
        addKey<VALUETYPE>(key);
        return values.insert( AnyPair_t(key, value) ).second;
      }

/**
 * Retrieve value associated with key
 *
 * @param key Key of value to return
 * @return value associated with key
 */
      template <typename VALUETYPE>
      VALUETYPE get(const std::string& key) const {
        AnyMap_t::const_iterator i = values.find(key);
        if (i != values.end())
          return Lucee::any_cast<VALUETYPE>((*i).second);
        // value not found: throw an exception
        Lucee::Except ex;
        ex << "Value " << key << " not found";
        throw ex;
      }

/**
 * Check if key exist in key-value pair.
 *
 * @param key Key of object to check.
 * @return true if key exists, false otherwise.
 */
      bool has(const std::string& key) const;

/**
 * Get list of keys for a particular type
 *
 * @return list of keys
 */
      template <typename VALUETYPE>
      std::vector<std::string> getKeys() const {
        return Lucee::typeMapExtract<VALUETYPE, TypeToKeys>(typeToKeys).keys;
      }

/** Typedef for map of strings to values, stored as Lucee::Any */
      typedef std::map<std::string, Lucee::Any, std::less<std::string> > AnyMap_t;
/** Typedef for pair of strings and values, stored as Lucee::Any */
      typedef std::pair<std::string, Lucee::Any> AnyPair_t;

    private:

/**
 * Container to maps types -> keys
 */
      template <typename T>
      struct TypeContainer {
/** List of keys of type T */
          std::vector<std::string> keys;
      };
/** Typedef for container for types to keys of those types */
      typedef Lucee::TypeMap<Lucee::DataTypes_t, TypeContainer> TypeToKeys;

/** Container mapping types to keys of those types */
      TypeToKeys typeToKeys;
/** Map of keys to values */
      AnyMap_t values;

/**
 * Add a key to type->keys map for supplied type
 *
 * @param key Key to add
 */
      template <typename VALUETYPE>
      void addKey(const std::string& key) {
        Lucee::typeMapExtract<VALUETYPE>(typeToKeys).keys.push_back(key);
      }

/**
 * Copy keys from give KeyVal object to this one
 *
 * @param kv KeyVal object to copy from
 */
      template <typename T>
      void copyKeys(const KeyVal& kv) {
        std::vector<std::string> nms = kv.getKeys<T>();
        std::vector<std::string>::const_iterator i;
        for (i=nms.begin(); i!=nms.end(); ++i)
          Lucee::typeMapExtract<T>(typeToKeys).keys.push_back(*i);
      }
  };

}

#endif // LC_KEY_VAL
