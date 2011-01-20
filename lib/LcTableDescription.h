/**
 * @file	LcTableDescription.h
 *
 * @brief	Description of a single value in a Lua table.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_TABLE_DESCRIPTION_H
#define LC_TABLE_DESCRIPTION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcExcept.h>
#include <LcLuaTable.h>
#include <LcValueDescription.h>
#include <LcVectorDescription.h>

// etc includes
#include <loki/HierarchyGenerators.h>

// std includes
#include <iostream>
#include <map>
#include <string>

namespace Lucee
{
// struct for checking/getting stuff from Lua table
  template <typename T> struct LuaTableFetcher;

// struct for checking/getting stuff from Lua table
  template <> struct LuaTableFetcher<int>
  {
/**
 * Check if key exists in Lua table.
 *
 * @param key Key to check.
 * @param tbl Lua table to check.
 * @return true if key exists, false otherwise.
 */
      static bool has(const std::string& key, Lucee::LuaTable& tbl) 
      {
        return tbl.hasNumber(key);
      }

/**
 * Get value for key in Lua table.
 *
 * @param key Key to check.
 * @param tbl Lua table to check.
 * @return value in table.
 */
      static int get(const std::string& key, Lucee::LuaTable& tbl) 
      {
        return (int) tbl.getNumber(key);
      }

/**
 * Return string description of type.
 *
 * @return string description of type.
 */
      static std::string typeString()
      {
        return "integer";
      }
  };

// struct for checking/getting stuff from Lua table
  template <> struct LuaTableFetcher<double>
  {
/**
 * Check if key exists in Lua table.
 *
 * @param key Key to check.
 * @param tbl Lua table to check.
 * @return true if key exists, false otherwise.
 */
      static bool has(const std::string& key, Lucee::LuaTable& tbl)
      {
        return tbl.hasNumber(key);
      }

/**
 * Get value for key in Lua table.
 *
 * @param key Key to check.
 * @param tbl Lua table to check.
 * @return value in table.
 */
      static double get(const std::string& key, Lucee::LuaTable& tbl) 
      {
        return tbl.getNumber(key);
      }

/**
 * Return string description of type.
 *
 * @return string description of type.
 */
      static std::string typeString()
      {
        return "number";
      }
  };

// struct for checking/getting stuff from Lua table
  template <> struct LuaTableFetcher<std::string>
  {
/**
 * Check if key exists in Lua table.
 *
 * @param key Key to check.
 * @param tbl Lua table to check.
 * @return true if key exists, false otherwise.
 */
      static bool has(const std::string& key, Lucee::LuaTable& tbl)
      {
        return tbl.hasString(key);
      }

/**
 * Get value for key in Lua table.
 *
 * @param key Key to check.
 * @param tbl Lua table to check.
 * @return value in table.
 */
      static std::string get(const std::string& key, Lucee::LuaTable& tbl) 
      {
        return tbl.getString(key);
      }

/**
 * Return string description of type.
 *
 * @return string description of type.
 */
      static std::string typeString()
      {
        return "string";
      }
  };

/**
 * This class describes a Lua table that can appear in a Lua
 * script. It mainly served as a validation mechanism for the Lua
 * table and also for autmatically generating documentation for each
 * table.
 */
  class TableDescription
  {
    public:
/**
 * Create a new table description object.
 *
 * @param nm Name of table.
 */
      TableDescription(const std::string& nm);

/**
 * Check and set values from supplied Lua table.
 *
 * @param tbl Lua table to set from.
 */
      void checkAndSet(Lucee::LuaTable& tbl);

/**
 * Add a new value to description. Reference to the value is returned
 * so that it can be specified in more detail if needed.
 *
 * @param nm Name of value to add.
 * @return reference to added value.
 */
      template <typename T>
      Lucee::ValueDescription<T>& addValue(const std::string& nm)
      {
        Lucee::ValueDescription<T> vd;
        Loki::Field<T>(vvTypeMap).values.insert(
          std::pair<std::string, Lucee::ValueDescription<T> >(nm, vd));
        typename std::map<std::string, Lucee::ValueDescription<T> >::iterator itr =
          Loki::Field<T>(vvTypeMap).values.find(nm);
        return itr->second;
      }

/**
 * Add a new optional value to description. Reference to the value is
 * returned so that it can be specified in more detail if needed.
 *
 * @param nm Name of value to add.
 * @param def Default value.
 * @return reference to added value.
 */
      template <typename T>
      Lucee::ValueDescription<T>& addValue(const std::string& nm, const T& def)
      {
        Lucee::ValueDescription<T> vd(def);
        Loki::Field<T>(vvTypeMap).values.insert(
          std::pair<std::string, Lucee::ValueDescription<T> >(nm, vd));
        typename std::map<std::string, Lucee::ValueDescription<T> >::iterator itr =
          Loki::Field<T>(vvTypeMap).values.find(nm);
        return itr->second;
      }

/**
 * Add a new value to description. Reference to the value is returned
 * so that it can be specified in more detail if needed.
 *
 * @param nm Name of value to add.
 * @return reference to added value.
 */
      template <typename T>
      Lucee::VectorDescription<T>& addVector(const std::string& nm)
      {
        Lucee::VectorDescription<T> vd;
        Loki::Field<T>(vvTypeMap).vectors.insert(
          std::pair<std::string, Lucee::VectorDescription<T> >(nm, vd));
        typename std::map<std::string, Lucee::VectorDescription<T> >::iterator itr =
          Loki::Field<T>(vvTypeMap).vectors.find(nm);
        return itr->second;
      }

/**
 * Add a new optional vector to description. Reference to the vector is
 * returned so that it can be specified in more detail if needed.
 *
 * @param nm Name of vector to add.
 * @param def Default vector.
 * @return reference to added vector.
 */
      template <typename T>
      Lucee::VectorDescription<T>& addVector(const std::string& nm, const T& def)
      {
        Lucee::VectorDescription<T> vd(def);
        Loki::Field<T>(vvTypeMap).vectors.insert(
          std::pair<std::string, Lucee::VectorDescription<T> >(nm, vd));
        typename std::map<std::string, Lucee::VectorDescription<T> >::iterator itr =
          Loki::Field<T>(vvTypeMap).vectors.find(nm);
        return itr->second;
      }

/**
 * Get value stored in table.
 *
 * @param nm Name of value to fetch.
 * @return value.
 */
      template <typename T>
      const Lucee::ValueDescription<T>& getValue(const std::string& nm) const
      {
        typename std::map<std::string, Lucee::ValueDescription<T> >::const_iterator itr =
          Loki::Field<T>(vvTypeMap).values.find(nm);
        if (itr == Loki::Field<T>(vvTypeMap).values.end())
        {
          Lucee::Except lce("TableDescription::getValue: Unable to find value '");
          lce << nm << "' in table" << std::endl;
          throw lce;
        }
        return itr->second;
      }

    private:
/** Name of table */
      std::string name;

/** Container for values and vectors of all types in table */
      template <typename T>
      struct ValVecContainer
      {
/** Value appearing in table */
          std::map<std::string, Lucee::ValueDescription<T> > values;
/** Vectors appearing in table */
          std::map<std::string, Lucee::VectorDescription<T> > vectors;
      };
/** Typedef to make life easier */
      typedef Loki::GenScatterHierarchy<LOKI_TYPELIST_3(int, double, std::string), ValVecContainer> ValVecTypeMap;
/** Container of values and vectors in table */
      ValVecTypeMap vvTypeMap;

/**
 * Check and set values from Lua table.
 *
 * @param tbl Table to set from.
 */
      template <typename T>
      bool checkAndSetValues(Lucee::LuaTable& tbl, std::ostringstream& errMsg)
      {
        bool pass=true;
// loop over each value
        typename std::map<std::string, Lucee::ValueDescription<T> >::iterator itr
          = Loki::Field<T>(vvTypeMap).values.begin();
        for (; itr != Loki::Field<T>(vvTypeMap).values.end(); ++itr)
        {
          if (LuaTableFetcher<T>::has(itr->first, tbl))
          { // value exists in table, check and fill it with table value
            T val = LuaTableFetcher<T>::get(itr->first, tbl);
            std::pair<bool, std::string> se = itr->second.checkValue(val);
            if (se.first)
              itr->second.fillVarWithValue(val);
            else
            { // check for value failed, report error
              pass = false;
              errMsg << "Validity test for '" << itr->first << "' failed." << std::endl;
              errMsg << "** ("  << se.second << ")" << std::endl;
            }
          }
          else
          { // value does not exist in table

            if (itr->second.isOptional())
            { // okay if missing and optional
              itr->second.fillVarWithOptional();
            }
            else
            { // error: required value not in table
              pass = false;
              errMsg << "Required " << LuaTableFetcher<T>::typeString()
                     << " '" << itr->first << "' not specifed" << std::endl;
            }
          }
        }
        return pass;
      }

/**
 * Check and set vector from Lua table.
 *
 * @param tbl Table to set from.
 */
      template <typename T>
      bool checkAndSetVectors(Lucee::LuaTable& tbl, std::ostringstream& errMsg)
      {
      }
  };
}

#endif // LC_TABLE_DESCRIPTION_H
