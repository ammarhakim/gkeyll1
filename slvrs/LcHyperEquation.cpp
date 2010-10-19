/**
 * @file	LcHyperEquation.cpp
 *
 * @brief	Interface to hyperbolic equations.
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
#include <LcExcept.h>
#include <LcHyperEquation.h>

namespace Lucee
{
// set module name
  const char *HyperEquation::id = "HyperEquation";

  HyperEquation::HyperEquation(unsigned meqn, unsigned mwave)
    : meqn(meqn), mwave(mwave)
  {
  }

  void
  HyperEquation::readInput(Lucee::LuaTable& tbl)
  {
  }

  void
  HyperEquation::flux(const Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double>& f)
  {
    throw Lucee::Except("HyperEquation::flux: Method not implemented");
  }
}
