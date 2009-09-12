/**
 * @file	LcKeyVal.h
 *
 * @brief	Class to hold key-value pairs.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

#ifndef LC_KEY_VAL_H
#define LC_KEY_VAL_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcExcept.h>

// Loki includes
#include <loki/HierarchyGenerators.h>

// std includes
#include <map>
#include <string>
#include <vector>

namespace Lucee
{
  typedef LOKI_TYPELIST_6(
      int,
      double,
      std::string,
      std::vector<int>,
      std::vector<double>,
      std::vector<std::string>) BasicTypeList_t;

/**
 * Lucee::KeyVal provides a container to store/retrive key-value
 * pairs. The key are strings, while the values can int, double,
 * std::string and vectors of these.
 *
 * Lucee::KeyVal objects are immutable, i.e, values once added can not
 * be removed or modifed.
 */
  class KeyVal
  {
    public:
/**
 * Create a new key-value object
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
 * Assignment operator: make a deep copy.
 *
 * @param kv KeyVal object to copy.
 * @return reference to assigned object.
 */
      KeyVal& operator=(const KeyVal& kv);

/**
 * Insert a key, value pair. Returns true if the insertion worked,
 * false otherwise.
 *
 * @param key Key of object.
 * @param value Value of object.
 * @return true, if insertion worked, false otherwise.
 */
      template<typename VALUETYPE>
      bool add(const std::string& key, VALUETYPE value) 
      {
        return Loki::Field<VALUETYPE>(dataStoreTypeMap).dataMap.insert(
          std::pair<std::string, VALUETYPE>(key, value)).second;
      }

/**
 * Insert a key, vector of values pair. Returns true if the insertion
 * worked, false otherwise.
 *
 * @param key Key of object.
 * @param value Value of object.
 * @return true, if insertion worked, false otherwise.
 */
      template<typename VALUETYPE>
      bool addVec(const std::string& key, const std::vector<VALUETYPE>& value) 
      {
        return Loki::Field<std::vector<VALUETYPE> >(dataStoreTypeMap).dataMap.insert(
          std::pair<std::string, std::vector<VALUETYPE> >(key, value)).second;
      }

/**
 * Retrieve value associated with key.
 *
 * @param key Key of value to return.
 * @return value associated with key.
 */
      template <typename VALUETYPE>
      VALUETYPE get(const std::string& key) const
      {
        typename std::map<std::string, VALUETYPE>::const_iterator itr = 
          Loki::Field<VALUETYPE>(dataStoreTypeMap).dataMap.find(key);
        if (itr != Loki::Field<VALUETYPE>(dataStoreTypeMap).dataMap.end())
          return itr->second;
        Lucee::Except lce("KeyVal::get: Key ");
        lce << key << " does not exists in KeyVal object" << std::endl;
        throw lce;
      }

/**
 * Retrieve vector value associated with key.
 *
 * @param key Key of value to return.
 * @return value associated with key.
 */
      template <typename VALUETYPE>
      std::vector<VALUETYPE> getVec(const std::string& key) const
      {
        typename std::map<std::string, std::vector<VALUETYPE> >::const_iterator itr = 
          Loki::Field<std::vector<VALUETYPE> >(dataStoreTypeMap).dataMap.find(key);
        if (itr != Loki::Field<std::vector<VALUETYPE> >(dataStoreTypeMap).dataMap.end())
          return itr->second;
        Lucee::Except lce("KeyVal::getVec: Key ");
        lce << key << " does not exists in KeyVal object" << std::endl;
        throw lce;
      }

/**
 * Check if key exist in key-value pair.
 *
 * @param key Key of object to check.
 * @return true if key exists, false otherwise.
 */
      template <typename VALUETYPE>
      bool has(const std::string& key) const
      {
        typename std::map<std::string, VALUETYPE>::const_iterator itr = 
          Loki::Field<VALUETYPE>(dataStoreTypeMap).dataMap.find(key);
        if (itr != Loki::Field<VALUETYPE>(dataStoreTypeMap).dataMap.end())
          return true;
        return false;
      }

/**
 * Check if key exist in key-vector-value pair.
 *
 * @param key Key of object to check.
 * @return true if key exists, false otherwise.
 */
      template <typename VALUETYPE>
      bool hasVec(const std::string& key) const
      {
        typename std::map<std::string, std::vector<VALUETYPE> >::const_iterator itr = 
          Loki::Field<std::vector<VALUETYPE> >(dataStoreTypeMap).dataMap.find(key);
        if (itr != Loki::Field<std::vector<VALUETYPE> >(dataStoreTypeMap).dataMap.end())
          return true;
        return false;
      }

    private:
/**
 * Structure to hold map of keys to values.
 */
      template <typename T>
      struct DataStore
      {
/** Map of keys to values */
          std::map<std::string, T> dataMap;
      };
/** Type definition for map of types to data containers */
      typedef Loki::GenScatterHierarchy<BasicTypeList_t, DataStore> DataStoreMap;
/** Map for containers for name -> data */
      DataStoreMap dataStoreTypeMap;
  };
}

#endif // LC_KEY_VAL_H
