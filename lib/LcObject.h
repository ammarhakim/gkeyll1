/**
 * @file	lcobject.h
 *
 * @brief	Base class for all Lucee objects
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

#ifndef LC_OBJECT_H
#define LC_OBJECT_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// std includes
#include <string>

namespace Lucee
{
/**
 * Base class for all Lucee objects.
 */
  class Object
  {
    public:
/**
 * Create a new object with given name.
 *
 * @param name Name of object.
 */
      Object(const std::string& name);

/**
 * Destroy object.
 */
      virtual ~Object();

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
      std::string name;
  };
}

#endif // LC_OBJECT_H
