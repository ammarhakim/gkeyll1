/**
 * @file	LcBasicTypeList.h
 *
 * @brief	Typelist for basic Lucee types.
 *
 * @version	$Id: LcBasicTypeList.h 195 2009-09-27 00:00:03Z a.hakim777 $
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

#ifndef LC_BASIC_TYPE_LIST_H
#define LC_BASIC_TYPE_LIST_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// Loki includes
#include <loki/HierarchyGenerators.h>

// std includes
#include <map>
#include <string>
#include <vector>

namespace Lucee
{
/** Basic supported types */
  typedef LOKI_TYPELIST_6(
    int,
    double,
    std::string,
    std::vector<int>,
    std::vector<double>,
    std::vector<std::string>) BasicTypeList_t;
}

#endif // LC_BASIC_TYPE_LIST_H
