/**
 * @file	LcUnstructConnectivity.cpp
 *
 * @brief	Class holding connectivity information.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcUnstructConnectivity.h>

namespace Lucee
{
  void UnstructConnectivity::reset(unsigned nd)
  {
// clear existing arrays
    indices.clear();
    offsets.clear();
// resize offset
    offsets.resize(nd+1);
  }
}
