/**
 * @file	LcPoissonBracketEquation.cpp
 *
 * @brief	Interface to poisson bracket equations.
 */
// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcExcept.h>
#include <LcGlobals.h>
#include <LcPoissonBracketEquation.h>
#include <LcPointerHolder.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  // set module name
  const char *PoissonBracketEquation::id = "PoissonBracketEquation";

  PoissonBracketEquation::PoissonBracketEquation()
  {
  }

  void
  PoissonBracketEquation::readInput(Lucee::LuaTable& tbl)
  {
  }

  void
  PoissonBracketEquation::computeAlphaAtQuadNodes(const Eigen::MatrixXd& gradHamiltonian, const Eigen::MatrixXd& interpMat,
    const int idx[], Eigen::MatrixXd& alpha)
  {
    throw Lucee::Except("PoissonBracketEquation::computeAlphaAtQuadNodes: Method not implemented");
  }
}
