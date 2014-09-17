/**
 * @file	LcPoissonBracketAdvectionEquation4D.cpp
 *
 * @brief	
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcPoissonBracketAdvectionEquation4D.h>

// std includes
#include <algorithm>

namespace Lucee
{
  const char *PoissonBracketAdvectionEquation4D::id = "AdvectionEquation4D";

  PoissonBracketAdvectionEquation4D::PoissonBracketAdvectionEquation4D()
    : Lucee::PoissonBracketEquation()
  {
  }

  void
  PoissonBracketAdvectionEquation4D::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::PoissonBracketEquation::readInput(tbl);
  }

  void
  PoissonBracketAdvectionEquation4D::computeAlphaAtQuadNodes(const Eigen::MatrixXd& hamiltonian, const Eigen::MatrixXd& interpMat,
          const int idx[], Eigen::MatrixXd& alpha)
  {
    alpha.setZero(alpha.rows(), alpha.cols());
    alpha.row(0).setOnes();
  }
}
