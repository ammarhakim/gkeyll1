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

    Eigen::MatrixXd massMatrix(nlocal, nlocal);
    
    nodalBasis->getMassMatrix(massMatrixLucee);
    
    // Get interpolation matrix, gaussian quadrature points, and weights
    Lucee::Matrix<double> interpMatrixLucee(numQuadNodes, nlocal);

    Lucee::Matrix<double> gaussOrdinatesLucee(numQuadNodes, 3);

    gaussWeights = std::vector<double>(numQuadNodes);

    // Allocate Eigen matrices
    interpMatrix = Eigen::MatrixXd(numQuadNodes, nlocal);
    gaussOrdinates = Eigen::MatrixXd(numQuadNodes, 3);

    // Get the interpolation matrix for the volume quadrature points
    nodalBasis->getGaussQuadData(interpMatrixLucee, gaussOrdinatesLucee, gaussWeights);

    copyLuceeToEigen(massMatrixLucee, massMatrix);
    copyLuceeToEigen(interpMatrixLucee, interpMatrix);
    copyLuceeToEigen(gaussOrdinatesLucee, gaussOrdinates);

    // Compute and store inverse of mass matrix
    massMatrixInv = massMatrix.inverse();

    // Compute and store transpose of interpolation matrix
    interpMatrixTranspose = interpMatrix.transpose();
  }

  Lucee::UpdaterStatus 
  ElectrostaticPhiUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<1>& grid
      = this->getGrid<Lucee::StructuredGridBase<1> >();

    // Electron and ion densities
    const Lucee::Field<1, double>& nElcIn = this->getInp<Lucee::Field<1, double> >(0);
    const Lucee::Field<1, double>& nIonIn = this->getInp<Lucee::Field<1, double> >(1);
    // Electron thermal velocity squared vt(x)^2
    const Lucee::Field<1, double>& vtSqIn = this->getInp<Lucee::Field<1, double> >(2);
    // Need to multiply by 1/(kPerpRhoSquared) in the Lua script
    Lucee::Field<1, double>& phiOut = this->getOut<Lucee::Field<1, double> >(0);

    int nlocal = nodalBasis->getNumNodes();

    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion();

    Lucee::ConstFieldPtr<double> nElcPtr = nElcIn.createConstPtr();
    Lucee::ConstFieldPtr<double> nIonPtr = nIonIn.createConstPtr();
    Lucee::ConstFieldPtr<double> vtSqPtr = vtSqIn.createConstPtr();
    Lucee::FieldPtr<double> phiPtr = phiOut.createPtr();

    phiPtr = 0.0;

    // Loop over all cells
    for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
    {
      // Set inputs
      nElcIn.setPtr(nElcPtr, ix);
      nIonIn.setPtr(nIonPtr, ix);
      vtSqIn.setPtr(vtSqPtr, ix);
      // Set outputs
      phiOut.setPtr(phiPtr, ix);

      Eigen::VectorXd nVec(nlocal);

      // Compute phi(x) at nodal points, then interpolate to gaussian quadrature points
      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
        nVec(componentIndex) = ELECTRON_MASS*vtSqPtr[componentIndex]*(1-nElcPtr[componentIndex]/nIonPtr[componentIndex])
          /(ELEMENTARY_CHARGE*kPerpTimesRho*kPerpTimesRho);
      Eigen::VectorXd phiAtQuadPoints = interpMatrix*nVec;
      
      // Compute projection of phi(x) onto basis functions
      for (int componentIndex = 0; componentIndex < phiAtQuadPoints.rows(); componentIndex++)
        phiAtQuadPoints(componentIndex) = gaussWeights[componentIndex]*phiAtQuadPoints(componentIndex);
      Eigen::VectorXd phiWeights = massMatrixInv*interpMatrixTranspose*phiAtQuadPoints;

      // Copy into output pointer
      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
        phiPtr[componentIndex] = phiWeights(componentIndex);
    }

    return Lucee::UpdaterStatus();
  }

  void
  ElectrostaticPhiUpdater::declareTypes()
  {
    // takes three inputs (n_e(x), n_i(x), vThermElecSq(x))
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
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
