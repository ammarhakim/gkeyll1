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

  QuadratureRule::QuadratureRule(unsigned numNodes)
    : numNodes(numNodes)
  {
  }

  QuadratureRule::~QuadratureRule()
  {
  }

  void
  QuadratureRule::readInput(Lucee::LuaTable& tbl)
  {
// get number of nodes
    if (tbl.hasNumber("numNodes"))
      numNodes = (unsigned) tbl.getNumber("numNodes");
    else
      Lucee::Except("QuadratureRule::readInput: Must provide 'numNodes' specifying number of nodes");
  }
}
