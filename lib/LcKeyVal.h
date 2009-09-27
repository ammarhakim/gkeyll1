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
#include <loki/Functor.h>
#include <loki/HierarchyGenerators.h>
#include <loki/SmartPtr.h>

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
 * be removed or modifed. Copying or assigning from another KeyVal
 * object creates shallow copies. Use the duplicate() method if a deep
 * copy is needed.
 */
  class KeyVal
  {
    public:
/**
 * Create a new key -> value set.
 */
      KeyVal();

/**
 * Copy ctor: make a shallow copy of supplied KeyVal class
 *
 * @param kv KeyVal object to copy.
 */
      KeyVal(const KeyVal& kv);

/**
 * Assignment operator: make a shallow copy.
 *
 * @param kv KeyVal object to copy.
 * @return reference to assigned object.
 */
      KeyVal& operator=(const KeyVal& kv);

/**
 * Get the number of values stored for the specified type.
 */
      template <typename VALUETYPE>
      unsigned getNum() const
      {
        return Loki::Field<VALUETYPE>(mapCntr->dataStoreTypeMap).dataMap.size();
      }

/**
 * Set to first element in key, values pair map.
 */
      template <typename VALUETYPE>
      void setToFirst() const
      {
        Loki::Field<VALUETYPE>(dataStoreItrTypeMap).itr =
          Loki::Field<VALUETYPE>(mapCntr->dataStoreTypeMap).dataMap.begin();
      }

/**
 * Return value of current key, value pair and increment the iterator
 * to the next entry.
 *
 * @param kvp key, value stored as a std::pair
 */
      template <typename VALUETYPE>
      std::pair<std::string, VALUETYPE> getAndBump() const
      {
        typename KeyVal::DataStoreItr<VALUETYPE>& dsi 
          = Loki::Field<VALUETYPE>(dataStoreItrTypeMap);
        std::pair<std::string, VALUETYPE> kvp(dsi.itr->first, dsi.itr->second);
        dsi.itr++;
        return kvp;
      }

/**
 * Insert a key, value pair. Returns true if the insertion worked,
 * false otherwise.
 *
 * @param key Key of object.
 * @param value Value of object.
 * @return true, if insertion worked, false otherwise.
 */
      template <typename VALUETYPE>
      bool add(const std::string& key, const VALUETYPE& value) 
      {
        return Loki::Field<VALUETYPE>(mapCntr->dataStoreTypeMap).dataMap.insert(
          std::pair<std::string, VALUETYPE>(key, value)).second;
      }

/**
 * Insert a new functor object. This should be a callable object:
 * i.e. a function pointer or an object with operator() overloaded.
 *
 * @param nm name of function.
 * @param func Functor object.
 * @return true, if addition worked, false otherwise.
 */
      bool addFunc(const std::string& nm,
        Loki::Functor<std::vector<double>, LOKI_TYPELIST_1(const std::vector<double>&)> func);

/**
 * Set internal iterator to function.
 *
 * @param nm Name of function.
 */
      void setToFunc(const std::string& nm);

/**
 * Evaulate the current function. This can be called only after the
 * setToFunc() method has been called.
 *
 * @param inp Parameters to pass to function.
 * @return value of function.
 */
      std::vector<double> evalCurrentFunc(const std::vector<double>& inp);

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
          Loki::Field<VALUETYPE>(mapCntr->dataStoreTypeMap).dataMap.find(key);
        if (itr != Loki::Field<VALUETYPE>(mapCntr->dataStoreTypeMap).dataMap.end())
          return itr->second;
        Lucee::Except lce("KeyVal::get: Key ");
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
          Loki::Field<VALUETYPE>(mapCntr->dataStoreTypeMap).dataMap.find(key);
        if (itr != Loki::Field<VALUETYPE>(mapCntr->dataStoreTypeMap).dataMap.end())
          return true;
        return false;
      }

/**
 * Create a duplicate from this object.
 *
 * @return duplicate of this object.
 */
      KeyVal duplicate() const;

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
/** Type definition for functor */
      typedef Loki::Functor<std::vector<double>, LOKI_TYPELIST_1(const std::vector<double>&)> Functor_t;

/** Structure to hold map */
      struct DataStoreMapCntr 
      { // this is a struct so we can reference count it
/** Map for containers for name -> data */
          DataStoreMap dataStoreTypeMap;
/** Map for key -> functor objects */
          std::map<std::string,  Functor_t> functorMap;
      };
/** Container for map */
      Loki::SmartPtr<DataStoreMapCntr> mapCntr;

/**
 * Iterator over key value pairs.
 */
      template <typename T>
      struct DataStoreItr
      {
/** Iterator into map */
          typename std::map<std::string, T>::const_iterator itr;
      };
/** Type definition for map of types to data containers iterators */
      typedef Loki::GenScatterHierarchy<BasicTypeList_t, DataStoreItr> DataStoreItrMap;
/** Map for containers for name -> data iterators */
      mutable DataStoreItrMap dataStoreItrTypeMap;
/** */
      mutable std::map<std::string, Functor_t>::iterator functorItr;

/**
 * Copy stuff over from supplied set.
 *
 * @param kv Set to copy from.
 */
      template <typename VALUETYPE>
      void copyToSet(KeyVal& dest) const
      {
// loop and copy from supplied set
        this->setToFirst<VALUETYPE>();
        for (unsigned i=0; i<this->getNum<VALUETYPE>(); ++i)
        {
          std::pair<std::string, VALUETYPE> p = this->getAndBump<VALUETYPE>();
          dest.add<VALUETYPE>(p.first, p.second);
        }
      }
  };
}

#endif // LC_KEY_VAL_H
