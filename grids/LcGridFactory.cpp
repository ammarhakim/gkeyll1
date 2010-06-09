/**
 * @file	LcGridFactory.cpp
 *
 * @brief	Base class for factories to create grids.
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
#include <LcGenericFactory.h>
#include <LcGridIfc.h>

namespace Lucee
{
// instantiate factory base class for grids and set module name
  template class Lucee::GenericFactory<Lucee::GridIfc>;
  template <> const char *Lucee::GenericFactory<Lucee::GridIfc>::id = "Grid";
}
