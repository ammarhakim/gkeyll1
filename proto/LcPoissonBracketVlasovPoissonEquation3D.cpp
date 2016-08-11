/**
 * @file	LcPoissonBracketVlasovPoissonEquation3D.cpp
 *
 * @brief	Poisson bracket for Vlasov-Poisson problem with 1D2V and static constant B.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcPoissonBracketVlasovPoissonEquation3D.h>

namespace Lucee
{
  const char *PoissonBracketVlasovPoissonEquation3D::id = "VlasovPoissonEquation3D";

  PoissonBracketVlasovPoissonEquation3D::PoissonBracketVlasovPoissonEquation3D()
    : Lucee::PoissonBracketEquation()
  {
  }

  void
  PoissonBracketVlasovPoissonEquation3D::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::PoissonBracketEquation::readInput(tbl);

    if (tbl.hasNumber("speciesMass"))
      speciesMass = tbl.getNumber("speciesMass");
    else
      throw Lucee::Except("LcPoissonBracketVlasovPoissonEquation3D: Must provide speciesMass.");

    if (tbl.hasNumber("speciesCharge"))
      speciesCharge = tbl.getNumber("speciesCharge");
    else
      throw Lucee::Except("LcPoissonBracketVlasovPoissonEquation3D: Must provide speciesCharge.");

    if (tbl.hasNumber("BZ0"))
      BZ0 = tbl.getNumber("BZ0");
    else
      throw Lucee::Except("LcPoissonBracketVlasovPoissonEquation3D: Must provide magnetic field constant BZ0.");
  }

  void
  PoissonBracketVlasovPoissonEquation3D::computeAlphaAtQuadNodes(const Eigen::MatrixXd& hamiltonian, const Eigen::MatrixXd& interpMat,
          const int idx[], Eigen::MatrixXd& alpha)
  {
    Eigen::Matrix3d poissonTensor = Eigen::Matrix3d::Zero();
    poissonTensor(0,1) = 1.0/speciesMass;
    poissonTensor(1,0) = -1.0/speciesMass;
    poissonTensor(1,2) = speciesCharge*BZ0/(speciesMass*speciesMass);
    poissonTensor(2,1) = -speciesCharge*BZ0/(speciesMass*speciesMass);
    alpha = poissonTensor*hamiltonian;
  }
}
