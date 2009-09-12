/**
 * @file	LcSettableObject.h
 *
 * @brief       Base class for all Lucee objects which can be set from name/values pairs.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

#ifndef LC_SETTABLE_OBJECT_H
#define LC_SETTABLE_OBJECT_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// Lucee includes
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
 * Base class for all Lucee objects that can be set from string/value
 * pairs.
 */
  class SettableObject
  {
    public:
/**
 * Create a new object with given name.
 *
 * @param name Name of object.
 */
      SettableObject(const std::string& name);

/**
 * Destroy object.
 */
      virtual ~SettableObject();

/**
 * Get object name.
 *
 * @return Name of object
 */
      std::string getName() const;

/**
 * Add a new data member that can be set.
 *
 * @param nm Name of data member to add.
 * @param data Pointer to the data object.
 * @param help Help message for this data object.
 */
      template <typename T>
      void addData(const std::string& nm, T *data, const std::string& help)
      {
// set data point and help message
        typename SettableData<T>::Data sd;
        sd.data = data;
        sd.help = help;
// insert it into the appropriate map
        Loki::Field<T>(settableTypeMap).dataMap[nm] = sd;
      }

/**
 * Set a new data member that can be set.
 *
 * @param nm Name of data member to set.
 * @param val Value to set.
 */
      template <typename T>
      void setData(const std::string& nm, const T& val)
      {
// fetch iterator to data object
        typename SettableData<T>::DataMap::iterator itr = 
          Loki::Field<T>(settableTypeMap).dataMap.find(nm);
        if (itr == Loki::Field<T>(settableTypeMap).dataMap.end())
        { // not found
          Lucee::Except lce("SettableData::setData: Data object ");
          lce << nm << std::endl;
          throw lce;
        }
        *(itr->second.data) = val;
      }

/**
 * Initialize object. This should set the object up so that it is in a
 * usable state. Should return true if the initialization was
 * successful, false otherwise.
 *
 * @return Success/failure status.
 */
      virtual bool init() = 0;

    private:
/**
 * Structure to hold help message and pointer to data.
 */
      template <typename T>
      struct SettableData
      {
          struct Data 
          {
/** Pointer to settable data */
              T *data;
/** Help message associated with data */
              std::string help;
          };
          typedef std::map<std::string, Data> DataMap;
          DataMap dataMap;
      };

/** Name of object */
      std::string name;
/** Type definition for map of types to I/O containers */
      typedef Loki::GenScatterHierarchy<BasicTypeList_t, SettableData> SettableTypeMap;
/** Map for containers for name -> data */
      SettableTypeMap settableTypeMap;
  };
}

#endif // LC_SETTABLE_OBJECT_H
