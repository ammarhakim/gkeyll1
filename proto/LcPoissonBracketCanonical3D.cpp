/**
 * @file	LcPoissonBracketCanonical3D.cpp
 *
 * @brief	
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcPoissonBracketCanonical3D.h>


namespace Lucee
{
  const char *PoissonBracketCanonical3D::id = "Canonical3D";

  PoissonBracketCanonical3D::PoissonBracketCanonical3D()
    : Lucee::PoissonBracketEquation()
  {
  }

  void
  PoissonBracketCanonical3D::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::PoissonBracketEquation::readInput(tbl);
  }

  void
  PoissonBracketCanonical3D::computeAlphaAtQuadNodes(const Eigen::MatrixXd& hamiltonian, const Eigen::MatrixXd& interpMat,
          const int idx[], Eigen::MatrixXd& alpha)
  {
    Eigen::Matrix3d poissonTensor = Eigen::Matrix3d::Zero();
    poissonTensor(0,1) = 1.0;
    poissonTensor(1,0) = -1.0;
    alpha = poissonTensor*hamiltonian;
  }
}
