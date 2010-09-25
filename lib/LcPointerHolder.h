/**
 * @file	LcPointerHolder.h
 *
 * @brief	Simple class to hold pointer
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_POINTER_HOLDER_H
#define LC_POINTER_HOLDER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

namespace Lucee
{
  template <typename T>
  struct PointerHolder
  {
      T *pointer;
  };

  template <typename T>
  void deletePtr(T *ptr)
  {
    delete ptr;
  }
}

#endif // LC_POINTER_HOLDER_H
