/**
 * @file	LcRteRegistry.cpp
 *
 * @brief	Class for registering RTE solver object.
 *
 * @version	$Id: LcRteRegistry.cpp 321 2010-03-04 23:13:16Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// lucee includes
#include <LcRteHomogeneousSlab.h>
#include <LcRteRegistry.h>
#include <LcSolverIfc.h>

namespace Lucee
{
  void registerRteObjects()
  {
// register stuff
    new Lucee::ObjRegistry<Lucee::SolverIfc, Lucee::RteHomogeneousSlab>;

  }
}

