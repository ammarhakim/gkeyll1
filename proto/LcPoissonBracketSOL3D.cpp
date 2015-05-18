/**
 * @file	LcPoissonBracketSOL3D.cpp
 *
 * @brief	Poisson bracket for SOL problem with 1D2V.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcPoissonBracketSOL3D.h>


namespace Lucee
{
  const char *PoissonBracketSOL3D::id = "SOL3D";

  PoissonBracketSOL3D::PoissonBracketSOL3D()
    : Lucee::PoissonBracketEquation()
  {
  }

  void
  PoissonBracketSOL3D::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::PoissonBracketEquation::readInput(tbl);

    if (tbl.hasNumber("speciesMass"))
      speciesMass = tbl.getNumber("speciesMass");
    else
      throw Lucee::Except("LcPoissonBracketSOL3D: Must provide speciesMass.");
  }

  void
  PoissonBracketSOL3D::computeAlphaAtQuadNodes(const Eigen::MatrixXd& hamiltonian, const Eigen::MatrixXd& interpMat,
          const int idx[], Eigen::MatrixXd& alpha)
  {
    Eigen::Matrix3d poissonTensor = Eigen::Matrix3d::Zero();
    poissonTensor(0,1) = 1/speciesMass;
    poissonTensor(1,0) = -1/speciesMass;
    alpha = poissonTensor*hamiltonian;
  }
}
