/**
 * @file	LcPoissonBracketAdvectionEquation5D.cpp
 *
 * @brief	
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcPoissonBracketAdvectionEquation5D.h>

// std includes
#include <algorithm>

namespace Lucee
{
  const char *PoissonBracketAdvectionEquation5D::id = "AdvectionEquation5D";

  PoissonBracketAdvectionEquation5D::PoissonBracketAdvectionEquation5D()
    : Lucee::PoissonBracketEquation()
  {
  }

  void
  PoissonBracketAdvectionEquation5D::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::PoissonBracketEquation::readInput(tbl);
  }

  void
  PoissonBracketAdvectionEquation5D::computeAlphaAtQuadNodes(const Eigen::MatrixXd& hamiltonian, const Eigen::MatrixXd& interpMat,
          const int idx[], Eigen::MatrixXd& alpha)
  {
    alpha = hamiltonian;
  }
}
