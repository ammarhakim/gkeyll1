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
 * Initialize object. This should set the object up so that it is in a
 * usable state. Should return true if the initialization was
 * successful, false otherwise.
 *
 * @return Success/failure status.
 */
      virtual bool init() = 0;

    private:
/** Name of object */
      std::string name;
/** Map of names -> integers */
      std::map<std::string, int*> intMap;
/** Map of names -> double */
      std::map<std::string, double*> doubleMap;
/** Map of names -> string */
      std::map<std::string, std::string*> strMap;
  };
}

#endif // LC_SETTABLE_OBJECT_H
