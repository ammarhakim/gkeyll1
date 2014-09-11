/**
 * @file	LcPoissonBracketGyroEquation4D.cpp
 *
 * @brief	
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcPoissonBracketGyroEquation4D.h>

// std includes
#include <algorithm>

namespace Lucee
{
  const char *PoissonBracketGyroEquation4D::id = "GyroEquation4D";

  PoissonBracketGyroEquation4D::PoissonBracketGyroEquation4D()
    : Lucee::PoissonBracketEquation()
  {

    //poissonTensor.resize(16);
  }

  void
  PoissonBracketGyroEquation4D::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::PoissonBracketEquation::readInput(tbl);

    if (tbl.hasObject<Lucee::Field<4, double> >("bParallelY"))
      bParallelYField = &tbl.getObject<Lucee::Field<4, double> >("bParallelY");
    else
      throw Lucee::Except("LcPoissonBracketGyroEquation4D: Must provide bParallelY.");
  }

  void
  PoissonBracketGyroEquation4D::computeAlphaAtQuadNodes(const Eigen::MatrixXd& hamiltonian, const Eigen::MatrixXd& interpMat,
          const int idx[], Eigen::MatrixXd& alpha)
  {
    int nlocal = interpMat.cols();
    Lucee::ConstFieldPtr<double> bParallelYPtr = bParallelYField->createConstPtr();
    bParallelYField->setPtr(bParallelYPtr, idx);
    
    double ionMass = 1.0;
    double c = 1.0;
    double elementaryCharge = 1.0;

    // Copy bParallelYField to an Eigen::VectorXd
    Eigen::VectorXd bParallelYVec(nlocal);
    for (int i = 0; i < nlocal; i++)
      bParallelYVec(i) = bParallelYPtr[i];
    
    alpha.setZero(alpha.rows(), alpha.cols());
    
    Eigen::Matrix4d poissonTensor = Eigen::Matrix4d::Zero();

    // (0, 1)
    Eigen::VectorXd poissonElement = -ionMass*ionMass*c/elementaryCharge*interpMat.rowwise().sum();

    alpha.row(0) += poissonElement.cwiseProduct(hamiltonian.row(1));

    // (1, 0)
    alpha.row(1) -= poissonElement.cwiseProduct(hamiltonian.row(0));

    // (1,2)
    poissonElement = ionMass*interpMat*bParallelYVec;
    alpha.row(1) += poissonElement.cwiseProduct(hamiltonian.row(2));

    // (2,1)
    alpha.row(2) -= poissonElement.cwiseProduct(hamiltonian.row(1));
  }
}
