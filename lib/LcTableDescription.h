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
#include <LcValueDescription.h>

// etc includes
#include <loki/HierarchyGenerators.h>

// std includes
#include <map>
#include <string>

namespace Lucee
{
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
 * Add a new value to description. Reference to the value is returned
 * so that it can be specified in more detail if needed.
 *
 * @param nm Name of value to add.
 * @return reference to added value.
 */
      template <typename T>
      Lucee::ValueDescription<T>& addValue(const std::string& nm)
      {
        Lucee::ValueDescription<T> vd(nm);
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
        Lucee::ValueDescription<T> vd(nm, def);
        Loki::Field<T>(vvTypeMap).values.insert(
          std::pair<std::string, Lucee::ValueDescription<T> >(nm, vd));
        typename std::map<std::string, Lucee::ValueDescription<T> >::iterator itr =
          Loki::Field<T>(vvTypeMap).values.find(nm);
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
      };
/** Typedef to make life easier */
      typedef Loki::GenScatterHierarchy<LOKI_TYPELIST_3(int, double, std::string), ValVecContainer> ValVecTypeMap;
/** Container of values and vectors in table */
      ValVecTypeMap vvTypeMap;
  };
}

#endif // LC_TABLE_DESCRIPTION_H
