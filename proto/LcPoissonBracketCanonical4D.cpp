/**
 * @file	LcPoissonBracketCanonical4D.cpp
 *
 * @brief	
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcPoissonBracketCanonical4D.h>


namespace Lucee
{
  const char *PoissonBracketCanonical4D::id = "Canonical4D";

  PoissonBracketCanonical4D::PoissonBracketCanonical4D()
    : Lucee::PoissonBracketEquation()
  {
  }

  void
  PoissonBracketCanonical4D::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::PoissonBracketEquation::readInput(tbl);
  }

  void
  PoissonBracketCanonical4D::computeAlphaAtQuadNodes(const Eigen::MatrixXd& hamiltonian, const Eigen::MatrixXd& interpMat,
          const int idx[], Eigen::MatrixXd& alpha)
  {
    Eigen::Matrix4d poissonTensor = Eigen::Matrix4d::Zero();
    //poissonTensor(0,2) = 1;
    //poissonTensor(1,3) = 1;
    //poissonTensor(2,0) = -1;
    //poissonTensor(3,1) = -1;
    poissonTensor(0,1) = 1;
    poissonTensor(1,0) = -1;
    alpha = poissonTensor*hamiltonian;
  }
}
