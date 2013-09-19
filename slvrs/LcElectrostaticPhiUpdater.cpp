/**
 * @file	LcElectrostaticPhiUpdater.cpp
 *
 * @brief	Updater to compute phi using a fixed value of k_perp*rho_s
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAlignedRectCoordSys.h>
#include <LcField.h>
#include <LcLinAlgebra.h>
#include <LcMathLib.h>
#include <LcElectrostaticPhiUpdater.h>
#include <LcStructuredGridBase.h>
#include <LcMathPhysConstants.h>
// for cutoff velocities
#include <LcDynVector.h>

// math include
#include <cmath>

namespace Lucee
{
// set id for module system
  const char *ElectrostaticPhiUpdater::id = "ElectrostaticPhiUpdater";

  ElectrostaticPhiUpdater::ElectrostaticPhiUpdater()
    : UpdaterIfc()
  {
  }

  void 
  ElectrostaticPhiUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis");
    else
      throw Lucee::Except("ElectrostaticPhiUpdater::readInput: Must specify element to use using 'basis'");

    if (tbl.hasNumber("kPerpTimesRho"))
      kPerpTimesRho = tbl.getNumber("kPerpTimesRho");
    else
      throw Lucee::Except("ElectrostaticPhiUpdater::readInput: Must specify kPerpTimesRho");

    if (tbl.hasNumber("Te0"))
      Te0 = tbl.getNumber("Te0");
    else
      throw Lucee::Except("ElectromagneticAUpdater::readInput: Must specify Te0");
  }

  void 
  ElectrostaticPhiUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();
    // global region to update
    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion();

    Lucee::RowMajorSequencer<1> seq(globalRgn);
    seq.step(); // just to get to first index
    int idx[1];
    seq.fillWithIndex(idx);
    nodalBasis->setIndex(idx);
    
    int nlocal = nodalBasis->getNumNodes();
    int numQuadNodes = nodalBasis->getNumGaussNodes();

    // Get mass matrix and then copy to Eigen format
    Lucee::Matrix<double> massMatrixLucee(nlocal, nlocal);

    massMatrix = Eigen::MatrixXd(nlocal, nlocal);
    
    nodalBasis->getMassMatrix(massMatrixLucee);
    
    // Get interpolation matrix, gaussian quadrature points, and weights
    Lucee::Matrix<double> interpMatrixLucee(numQuadNodes, nlocal);

    Lucee::Matrix<double> gaussOrdinatesLucee(numQuadNodes, 3);

    std::vector<double> gaussWeights(numQuadNodes);

    // Allocate Eigen matrices
    Eigen::MatrixXd interpMatrix(numQuadNodes, nlocal);
    Eigen::MatrixXd gaussOrdinates(numQuadNodes, 3);

    // Get the interpolation matrix for the volume quadrature points
    nodalBasis->getGaussQuadData(interpMatrixLucee, gaussOrdinatesLucee, gaussWeights);

    copyLuceeToEigen(massMatrixLucee, massMatrix);
    copyLuceeToEigen(interpMatrixLucee, interpMatrix);
    copyLuceeToEigen(gaussOrdinatesLucee, gaussOrdinates);

    tripleProducts = std::vector<Eigen::MatrixXd>(nlocal);
    // Create and store the triple-product basis integrals

    for (int basisIndex = 0; basisIndex < nlocal; basisIndex++)
    {
      // Initialize matrices
      tripleProducts[basisIndex]  = Eigen::MatrixXd::Zero(nlocal, nlocal);
      for (int rowIndex = 0; rowIndex < nlocal; rowIndex++)
      {
        for (int colIndex = 0; colIndex < nlocal; colIndex++)
        {
          double integralResult = 0.0;

          for (int gaussNodeIndex = 0; gaussNodeIndex < numQuadNodes; gaussNodeIndex++)
          {
            integralResult += gaussWeights[gaussNodeIndex]*interpMatrix(gaussNodeIndex, basisIndex)*
              interpMatrix(gaussNodeIndex, rowIndex)*interpMatrix(gaussNodeIndex, colIndex);
          }

          tripleProducts[basisIndex](rowIndex, colIndex) = integralResult;
        }
      }
    }
  }

  Lucee::UpdaterStatus 
  ElectrostaticPhiUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<1>& grid
      = this->getGrid<Lucee::StructuredGridBase<1> >();

    // Electron and ion densities
    const Lucee::Field<1, double>& nElcIn = this->getInp<Lucee::Field<1, double> >(0);
    const Lucee::Field<1, double>& nIonIn = this->getInp<Lucee::Field<1, double> >(1);
    // Dynvector containing the cutoff velocities computed at one or both edges
    const Lucee::DynVector<double>& cutoffVIn = this->getInp<Lucee::DynVector<double> >(2);
    Lucee::Field<1, double>& phiOut = this->getOut<Lucee::Field<1, double> >(0);

    int nlocal = nodalBasis->getNumNodes();

    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion();

    Lucee::ConstFieldPtr<double> nElcPtr = nElcIn.createConstPtr();
    Lucee::ConstFieldPtr<double> nIonPtr = nIonIn.createConstPtr();
    Lucee::FieldPtr<double> phiPtr = phiOut.createPtr();

    // Compute cutoff velocity on the right edge
    std::vector<double> cutoffVelocities = cutoffVIn.getLastInsertedData();
    double phiS = 0.5*ELECTRON_MASS*cutoffVelocities[1]*cutoffVelocities[1]/ELEMENTARY_CHARGE;

    phiPtr = 0.0;
    
    // Loop over all cells
    for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
    {
      // Set inputs
      nElcIn.setPtr(nElcPtr, ix);
      nIonIn.setPtr(nIonPtr, ix);
      // Set outputs
      phiOut.setPtr(phiPtr, ix);

      Eigen::VectorXd nIonVec(nlocal);
      Eigen::VectorXd nDiffVec(nlocal);
      
      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
      {
        nIonVec(componentIndex) = nIonPtr[componentIndex];
        nDiffVec(componentIndex) = nIonPtr[componentIndex] - nElcPtr[componentIndex];
        // Scale by various factors in the phi-equation
        nDiffVec(componentIndex) = Te0*ELEMENTARY_CHARGE*nDiffVec(componentIndex)/(ELEMENTARY_CHARGE*kPerpTimesRho*kPerpTimesRho);
      }

      Eigen::VectorXd rhsIntegrals = massMatrix*nDiffVec;

      Eigen::MatrixXd phiProjections(nlocal, nlocal);

      for (int basisIndex = 0; basisIndex < nlocal; basisIndex++)
      {
        Eigen::VectorXd phiProjectionsSingle = tripleProducts[basisIndex]*nIonVec;
        // Store components into rows of phiProjections matrix
        for (int colIndex = 0; colIndex < nlocal; colIndex++)
          phiProjections(basisIndex, colIndex) = phiProjectionsSingle(colIndex);
      }

      // Compute phi(x) weights
      Eigen::VectorXd phiWeights = phiProjections.inverse()*rhsIntegrals;

      // Copy into output pointer
      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
        phiPtr[componentIndex] = phiWeights(componentIndex) + phiS;
    }

    return Lucee::UpdaterStatus();
  }

  void
  ElectrostaticPhiUpdater::declareTypes()
  {
    // takes three inputs (n_e(x), n_i(x), vThermElecSq(x)) + cutoff velocities at edges
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::DynVector<double>));
    // returns one output: phi(x)
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }

  void
  ElectrostaticPhiUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}
