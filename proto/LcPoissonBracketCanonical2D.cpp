/**
 * @file	LcPoissonBracketCanonical2D.cpp
 *
 * @brief	
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcPoissonBracketCanonical2D.h>

// std includes
#include <algorithm>

namespace Lucee
{
  const char *PoissonBracketCanonical2D::id = "Canonical2D";

  PoissonBracketCanonical2D::PoissonBracketCanonical2D()
    : Lucee::PoissonBracketEquation()
  {
  }

  void
  PoissonBracketCanonical2D::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::PoissonBracketEquation::readInput(tbl);
  }

  void
  PoissonBracketCanonical2D::computeAlphaAtQuadNodes(const Eigen::MatrixXd& hamiltonian, const Eigen::MatrixXd& interpMat,
          const int idx[], Eigen::MatrixXd& alpha)
  {
    //alpha.setZero(alpha.rows(), alpha.cols());
    //alpha.row(0).setOnes();
    Eigen::Matrix2d poissonTensor = Eigen::Matrix2d::Zero();
    poissonTensor(0,1) = 1;
    poissonTensor(1,0) = -1;
    alpha = poissonTensor*hamiltonian;
  }
}
