/**
 * @file	lcdatatypes.h
 *
 * @brief	Class of supported types in I/O and messaging classes.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_DATA_TYPES_H
#define LC_DATA_TYPES_H

// lib includes
#include <lcany.h>

// etc includes
#include <loki/TypelistMacros.h>

// std includes
#include <string>
#include <vector>

namespace Lucee
{
/**
 * Typelist for supported I/O and message-ing types. These can be
 * augment if needed.
 *
 * NOTE: If adding more types change typelist length. Also ensure no
 * duplicates exist in the list.
 */
  typedef LOKI_TYPELIST_19(
      bool,
      char,
      unsigned char,
      short,
      unsigned short,
      int,
      unsigned int,
      long,
      unsigned long,
      float,
      double,
      long double,
      long long int,
      std::string, 
      std::vector<int>,
      std::vector<float>,
      std::vector<double>,
      std::vector<std::string>,
      std::vector<Lucee::Any>
                           ) DataTypes_t;
}

#endif // LC_DATA_TYPES_H
