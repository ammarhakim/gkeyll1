/**
 * @file	lckeyval.h
 *
 * @brief	Class to hold key-value pairs.
 *
 * @version	$Id$ *
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
 * Lucee::KeyVal objects are immutable: values once added can not be
 * removed or modifed.
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
      KeyVal& operator=(const KeyVal& rhs);

/**
 * Insert a name, value pair. Returns true if the insertion worked,
 * false otherwise
 *
 * @param name Name of object.
 * @param value Value of object.
 * @return true, if insertion worked, false otherwise.
 */
      template<typename VALUETYPE>
      bool add(const std::string& name, VALUETYPE value) {
        addName<VALUETYPE>(name);
        return _values.insert( AnyPair_t(name, value) ).second;
      }

/**
 * Check if key exist in key-value pair.
 *
 * @param key Key of object to check.
 * @return true if name exists, false otherwise.
 */
      bool has(const std::string& key) const;

/**
 * Retrieve value associated with key
 *
 * @param key Key of value to return
 * @return value associated with key
 */
      template <typename VALUETYPE>
      VALUETYPE get(const std::string& key) const {
        AnyMap_t::const_iterator i = _values.find(name);
        if (i != _values.end())
          return wx_any_cast<VALUETYPE>((*i).second);
        // value not found: throw an exception
        WxExcept wxe;
        wxe << "Value " << name << " not found";
        throw wxe;
      }

/**
 * Retrieve list of values associated with name
 *
 * @param name Return list of values associated with 'name'
 */
      template <typename VALUETYPE>
      std::vector<VALUETYPE> getVec(const std::string& name) const {
        // first fetch vector of WxAnys
        std::vector<WxAny> vals = 
          this->template get<std::vector<WxAny> >(name);
        // now convert it into vector of VALUETYPE
        std::vector<VALUETYPE> res;
        for (unsigned i=0; i<vals.size(); ++i)
          res.push_back( wx_any_cast<VALUETYPE>(vals[i]) );
        return res;
      }

/**
 * Get names of inserted types 
 *
 * @return list of names
 */
      template <typename VALUETYPE>
      std::vector<std::string> getNames() const {
        return wxTypeMapExtract<VALUETYPE, WxTypeToNames>(_typeToNames).names;
      }

/** Typedef for map of strings to values, stored as Lucee::Any */
      typedef std::map<std::string, Lucee::Any, std::less<std::string> > AnyMap_t;
/** Typedef for pair of strings and values, stored as Lucee::Any */
      typedef std::pair<std::string, Lucee::Any> AnyPair_t;

    private:
/**
 * Container to maps types -> names
 */
      template <typename T>
      struct TypeContainer {
          std::vector<std::string> names;
      };
      typedef WxTypeMap<WxDataTypes_t, TypeContainer> TypeToNames;

      TypeToNames _typeToNames;
      AnyMap_t _values;

/**
 * Add a name to type->names map for supplied type
 */
      template <typename VALUETYPE>
      void addName(const std::string& name) {
        wxTypeMapExtract<VALUETYPE>(_typeToNames).names.push_back(name);
      }

/**
 * Copy names from crypt type map to this object
 */
      template <typename T>
      void copyNames(const KeyVal& crypt) {
        std::vector<std::string> nms = crypt.getNames<T>();
        std::vector<std::string>::const_iterator i;
        for (i=nms.begin(); i!=nms.end(); ++i)
          wxTypeMapExtract<T>(_typeToNames).names.push_back(*i);
      }
  };

}

#endif // LC_KEY_VAL
