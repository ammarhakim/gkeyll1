/**
 * @file	lcdatatypes.h
 *
 * @brief	Class of supported types in I/O and messaging classes.
 *
 * @version	$Id$ *
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_DATA_TYPES_H
#define LC_DATA_TYPES_H

// lin includes
#include <lctypelist.h>
#include <lcany.h>

// std includes
#include <string>
#include <vector>

namespace Lucee
{
// Typelist for supported I/O and message-ing types. These can be
// augment if needed.
//
// NOTE: If adding more types change typelist length. Also ensure no
// duplicates exist in the list.
  typedef LC_TYPELIST_18(
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
      Lucee::Any,
      std::vector<Lucee::Any>,
      std::string,
      float*,
      double*
                         ) DataTypes_t;
}

#endif // LC_DATA_TYPES_H
