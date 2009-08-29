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

// std includes
#include <map>
#include <string>
#include <vector>

namespace Lucee
{
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
      void addData(const std::string& nm, int *data, const std::string& help);

/**
 * Add a new data member that can be set.
 *
 * @param nm Name of data member to add.
 * @param data Pointer to the data object.
 * @param help Help message for this data object.
 */
      void addData(const std::string& nm, double *data, const std::string& help);

/**
 * Add a new data member that can be set.
 *
 * @param nm Name of data member to add.
 * @param data Pointer to the data object.
 * @param help Help message for this data object.
 */
      void addData(const std::string& nm, std::string *data, const std::string& help);

/**
 * Set a new data member that can be set.
 *
 * @param nm Name of data member to set.
 * @param val Value to set.
 */
      void setData(const std::string& nm, int val);

/**
 * Set a new data member that can be set.
 *
 * @param nm Name of data member to set.
 * @param val Value to set.
 */
      void setData(const std::string& nm, double val);

/**
 * Set a new data member that can be set.
 *
 * @param nm Name of data member to set.
 * @param val Value to set.
 */
      void setData(const std::string& nm, const std::string& val);

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
/** Pointer to settable data */
          T *data;
/** Help message associated with data */
          std::string help;
      };

/** Name of object */
      std::string name;
/** Map of names -> integers */
      std::map<std::string, SettableData<int> > intMap;
/** Map of names -> double */
      std::map<std::string, SettableData<double> > doubleMap;
/** Map of names -> string */
      std::map<std::string, SettableData<std::string> > strMap;
  };
}

#endif // LC_SETTABLE_OBJECT_H
