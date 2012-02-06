/**
 * @file	LcQuadratureRule.cpp
 *
 * @brief	Base class for a quadrature rule.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcQuadratureRule.h>

namespace Lucee
{
// set module name
  const char *QuadratureRule::id = "QuadratureRule";

  void
  QuadratureRule::readInput(Lucee::LuaTable& tbl)
  {
  }

  QuadratureRule::~QuadratureRule()
  {
  }
}
