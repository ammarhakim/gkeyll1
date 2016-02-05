/**
 * @file	LcPoissonBracketGyroEquation5D.cpp
 *
 * @brief	Compute Jacobian*alpha at all quadrature points
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcPoissonBracketGyroEquation5D.h>

namespace Lucee
{
  const char *PoissonBracketGyroEquation5D::id = "GyroEquation5D";

  PoissonBracketGyroEquation5D::PoissonBracketGyroEquation5D()
    : Lucee::PoissonBracketEquation()
  {
  }

  void
  PoissonBracketGyroEquation5D::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::PoissonBracketEquation::readInput(tbl);

    if (tbl.hasObject<Lucee::Field<5, double> >("bStarY"))
      bStarYField = &tbl.getObject<Lucee::Field<5, double> >("bStarY");
    else
      throw Lucee::Except("LcPoissonBracketGyroEquation5D: Must provide bStarY.");

    if (tbl.hasObject<Lucee::Field<5, double> >("bStarZ"))
      bStarZField = &tbl.getObject<Lucee::Field<5, double> >("bStarZ");
    else
      throw Lucee::Except("LcPoissonBracketGyroEquation5D: Must provide bStarZ.");

    if (tbl.hasNumber("speciesCharge"))
      speciesCharge = tbl.getNumber("speciesCharge");
    else
      throw Lucee::Except("LcPoissonBracketGyroEquation5D: Must provide speciesCharge.");

    if (tbl.hasNumber("speciesMass"))
      speciesMass = tbl.getNumber("speciesMass");
    else
      throw Lucee::Except("LcPoissonBracketGyroEquation5D: Must provide speciesMass.");
  }

  void
  PoissonBracketGyroEquation5D::computeAlphaAtQuadNodes(const Eigen::MatrixXd& hamiltonian, const Eigen::MatrixXd& interpMat,
          const int idx[], Eigen::MatrixXd& alpha)
  {
    int nlocal = interpMat.cols();

    Lucee::ConstFieldPtr<double> bStarYPtr = bStarYField->createConstPtr();
    Lucee::ConstFieldPtr<double> bStarZPtr = bStarZField->createConstPtr();
    
    bStarYField->setPtr(bStarYPtr, idx);
    bStarZField->setPtr(bStarZPtr, idx);
    // Copy bStarYField and bStarZField to an Eigen::VectorXd
    Eigen::VectorXd bStarYVec(nlocal);
    Eigen::VectorXd bStarZVec(nlocal);

    for (int i = 0; i < nlocal; i++)
    {
      bStarYVec(i) = bStarYPtr[i];
      bStarZVec(i) = bStarZPtr[i];
    }
    alpha.setZero(alpha.rows(), alpha.cols());

    Eigen::RowVectorXd poissonElement = -speciesMass*speciesMass/speciesCharge*interpMat.rowwise().sum();
    // (0, 1)
    alpha.row(0) = poissonElement.cwiseProduct(hamiltonian.row(1));
    // (1, 0)
    alpha.row(1) = -poissonElement.cwiseProduct(hamiltonian.row(0));
    
    // Get a vector of bStarYVec*speciesMass at all quadrature points
    poissonElement = speciesMass*interpMat*bStarYVec;
    // (1,3)
    alpha.row(1) += poissonElement.cwiseProduct(hamiltonian.row(3));
    // (3,1)
    alpha.row(3) = -poissonElement.cwiseProduct(hamiltonian.row(1));

    // Get a vector of bStarZVec*speciesMass at all quadrature points
    poissonElement = speciesMass*interpMat*bStarZVec;
    // (2,3)
    alpha.row(2) = poissonElement.cwiseProduct(hamiltonian.row(3));
    // (3,2)
    alpha.row(3) -= poissonElement.cwiseProduct(hamiltonian.row(2));
  }

  void
  PoissonBracketGyroEquation5D::computeAlphaAtQuadNodes(const Eigen::MatrixXd& hamiltonian, const Eigen::MatrixXd& interpMat,
    const std::vector<int>& nodeNums, const int idx[], Eigen::RowVectorXd& alphaDotN, const int component)
  {
    int nlocal = nodeNums.size();
    //interpMat.cols();

    Lucee::ConstFieldPtr<double> bStarYPtr = bStarYField->createConstPtr();
    Lucee::ConstFieldPtr<double> bStarZPtr = bStarZField->createConstPtr();
    
    bStarYField->setPtr(bStarYPtr, idx);
    bStarZField->setPtr(bStarZPtr, idx);
    // Copy bStarYField and bStarZField to an Eigen::VectorXd
    Eigen::VectorXd bStarYVec(nlocal);
    Eigen::VectorXd bStarZVec(nlocal);

    for (int i = 0; i < nlocal; i++)
    {
      bStarYVec(i) = bStarYPtr[nodeNums[i]];
      bStarZVec(i) = bStarZPtr[nodeNums[i]];
    }

    Eigen::RowVectorXd poissonElement;

    // Only return a single row of the complete alpha matrix (gives a single component of alpha
    // evaluated at requested quadrature points)
    if (component == 0)
    {
      poissonElement = -speciesMass*speciesMass/speciesCharge*interpMat.rowwise().sum();
      alphaDotN = poissonElement.cwiseProduct(hamiltonian.row(1));
    }
    else if (component == 1)
    {
      poissonElement = -speciesMass*speciesMass/speciesCharge*interpMat.rowwise().sum();
      alphaDotN = -poissonElement.cwiseProduct(hamiltonian.row(0));
      poissonElement = speciesMass*interpMat*bStarYVec;
      alphaDotN += poissonElement.cwiseProduct(hamiltonian.row(3));
    }
    else if (component == 2)
    {
      poissonElement = speciesMass*interpMat*bStarZVec;
      alphaDotN = poissonElement.cwiseProduct(hamiltonian.row(3));
    }
    else if (component == 3)
    {
      poissonElement = speciesMass*interpMat*bStarYVec;
      alphaDotN = -poissonElement.cwiseProduct(hamiltonian.row(1));
      poissonElement = speciesMass*interpMat*bStarZVec;
      alphaDotN -= poissonElement.cwiseProduct(hamiltonian.row(2));
    }
  }
}
