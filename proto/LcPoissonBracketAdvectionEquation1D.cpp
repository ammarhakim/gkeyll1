/**
 * @file	LcPoissonBracketAdvectionEquation1D.cpp
 *
 * @brief	
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcPoissonBracketAdvectionEquation1D.h>

// std includes
#include <algorithm>

namespace Lucee
{
  const char *PoissonBracketAdvectionEquation1D::id = "AdvectionEquation1D";

  PoissonBracketAdvectionEquation1D::PoissonBracketAdvectionEquation1D()
    : Lucee::PoissonBracketEquation()
  {
  }

  void
  PoissonBracketAdvectionEquation1D::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::PoissonBracketEquation::readInput(tbl);
  }

  void
  PoissonBracketAdvectionEquation1D::computeAlphaAtQuadNodes(const Eigen::MatrixXd& hamiltonian, const Eigen::MatrixXd& interpMat,
          const int idx[], Eigen::MatrixXd& alpha)
  {
    //alpha.setZero(alpha.rows(), alpha.cols());
    //alpha.row(0).setOnes();
    alpha = hamiltonian;
  }
}
