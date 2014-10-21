/**
 * @file	LcPoissonBracketGyroEquation4D.cpp
 *
 * @brief	Compute Jacobian*alpha at all quadrature points
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcPoissonBracketGyroEquation4D.h>

namespace Lucee
{
  const char *PoissonBracketGyroEquation4D::id = "GyroEquation4D";

  PoissonBracketGyroEquation4D::PoissonBracketGyroEquation4D()
    : Lucee::PoissonBracketEquation()
  {
  }

  void
  PoissonBracketGyroEquation4D::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::PoissonBracketEquation::readInput(tbl);

    if (tbl.hasObject<Lucee::Field<4, double> >("bStarY"))
      bStarYField = &tbl.getObject<Lucee::Field<4, double> >("bStarY");
    else
      throw Lucee::Except("LcPoissonBracketGyroEquation4D: Must provide bStarY.");

    if (tbl.hasNumber("speciesCharge"))
      speciesCharge = tbl.getNumber("speciesCharge");
    else
      throw Lucee::Except("LcPoissonBracketGyroEquation4D: Must provide speciesCharge.");

    if (tbl.hasNumber("speciesMass"))
      speciesMass = tbl.getNumber("speciesMass");
    else
      throw Lucee::Except("LcPoissonBracketGyroEquation4D: Must provide speciesMass.");
  }

  void
  PoissonBracketGyroEquation4D::computeAlphaAtQuadNodes(const Eigen::MatrixXd& hamiltonian, const Eigen::MatrixXd& interpMat,
          const int idx[], Eigen::MatrixXd& alpha)
  {
    int nlocal = interpMat.cols();
    Lucee::ConstFieldPtr<double> bStarYPtr = bStarYField->createConstPtr();
    
    bStarYField->setPtr(bStarYPtr, idx);
    // Copy bStarYField to an Eigen::VectorXd
    Eigen::VectorXd bStarYVec(nlocal);
    for (int i = 0; i < nlocal; i++)
      bStarYVec(i) = bStarYPtr[i];

    alpha.setZero(alpha.rows(), alpha.cols());

    Eigen::RowVectorXd poissonElement = -speciesMass*speciesMass/speciesCharge*interpMat.rowwise().sum();
    // (0, 1)
    alpha.row(0) += poissonElement.cwiseProduct(hamiltonian.row(1));
    // (1, 0)
    alpha.row(1) -= poissonElement.cwiseProduct(hamiltonian.row(0));
    
    // Get a vector of bStarYVec*speciesMass at all quadrature points
    poissonElement = speciesMass*interpMat*bStarYVec;
    // (1,2)
    alpha.row(1) += poissonElement.cwiseProduct(hamiltonian.row(2));
    // (2,1)
    alpha.row(2) -= poissonElement.cwiseProduct(hamiltonian.row(1));
  }
}
