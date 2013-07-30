/**
 * @file	LcBoltzmannPhiUpdater.cpp
 *
 * @brief	Updater to compute average v_therm(x)^2 times ln(ni(x)/ni0)
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
#include <LcBoltzmannPhiUpdater.h>
#include <LcStructuredGridBase.h>

// math include
#include <cmath>

namespace Lucee
{
// set id for module system
  const char *BoltzmannPhiUpdater::id = "BoltzmannPhiUpdater";

  BoltzmannPhiUpdater::BoltzmannPhiUpdater()
    : UpdaterIfc()
  {
  }

  void 
  BoltzmannPhiUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasNumber("electronMass"))
      electronMass = tbl.getNumber("electronMass");
    else
      throw Lucee::Except("BoltzmannPhiUpdater::readInput: Must specify electronMass");

    if (tbl.hasNumber("ionMass"))
      ionMass = tbl.getNumber("ionMass");
    else
      throw Lucee::Except("BoltzmannPhiUpdater::readInput: Must specify ionMass");

    if (tbl.hasNumber("elementaryCharge"))
      elementaryCharge = tbl.getNumber("elementaryCharge");
    else
      throw Lucee::Except("BoltzmannPhiUpdater::readInput: Must specify elementaryCharge");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis");
    else
      throw Lucee::Except("BoltzmannPhiUpdater::readInput: Must specify element to use using 'basis'");
  }

  void 
  BoltzmannPhiUpdater::initialize()
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
  BoltzmannPhiUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<1>& grid
      = this->getGrid<Lucee::StructuredGridBase<1> >();

    // Particle density n(x)
    const Lucee::Field<1, double>& nIn = this->getInp<Lucee::Field<1, double> >(0);
    // First velocity moment of ion distribution function
    const Lucee::Field<1, double>& mom1In = this->getInp<Lucee::Field<1, double> >(1);
    // Thermal velocity squared vt(x)^2
    const Lucee::Field<1, double>& vtSqIn = this->getInp<Lucee::Field<1, double> >(2);
    // Returns (e/m_i)*phi(x), so need to multiply by m_i/e outside of code
    Lucee::Field<1, double>& phiOut = this->getOut<Lucee::Field<1, double> >(0);

    int nlocal = nodalBasis->getNumNodes();

    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion();

    Lucee::ConstFieldPtr<double> nPtr = nIn.createConstPtr();
    Lucee::ConstFieldPtr<double> mom1Ptr = mom1In.createConstPtr();
    Lucee::ConstFieldPtr<double> vtSqPtr = vtSqIn.createConstPtr();
    Lucee::FieldPtr<double> phiPtr = phiOut.createPtr();

    phiPtr = 0.0;
    
    // Find value of n(x) at the very last edge of the domain
    nIn.setPtr(nPtr, globalRgn.getUpper(0)-1);
    mom1In.setPtr(mom1Ptr, globalRgn.getUpper(0)-1);
    double nEdge = nPtr[nlocal-1];
    double ionFluxAtEdge = mom1Ptr[nlocal-1];

    // Find mean vtSq in entire domain
    double meanVtSq = 0.0;

    for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
    {
      // Set inputs
      vtSqIn.setPtr(vtSqPtr, ix);
      Eigen::VectorXd vtSqVec(nlocal);
      
      // Figure out vtSq(x) at quadrature points in the cell
      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
        vtSqVec(componentIndex) = vtSqPtr[componentIndex];
      Eigen::VectorXd vtSqAtQuadPoints = interpMatrix*vtSqVec;

      for (int componentIndex = 0; componentIndex < vtSqAtQuadPoints.rows(); componentIndex++)
        meanVtSq += gaussWeights[componentIndex]*vtSqAtQuadPoints(componentIndex);
    }

    // Divide by length of domain
    // Consider using grid.getNumCells(0)
    meanVtSq = meanVtSq/(grid.getDx(0)*(globalRgn.getUpper(0)-globalRgn.getLower(0)));

    // Electron temperature in joules
    double Te = ionMass*meanVtSq;
    // Electron thermal velocity
    double vThermElec = sqrt(Te/electronMass);
    const double PI = 3.141592653589793238462;
    // Figure out phi_s (sheath potential)
    double phiS = -Te/elementaryCharge*log(sqrt(2*PI)*ionFluxAtEdge/(nEdge*vThermElec));

    if (nEdge < 0.0)
    {
      throw Lucee::Except("BoltzmannPhiUpdater::update: nEdge is negative!");
    }


    // Loop over all cells
    for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
    {
      // Set inputs
      nIn.setPtr(nPtr, ix);
      // Set outputs
      phiOut.setPtr(phiPtr, ix);

      Eigen::VectorXd nVec(nlocal);

      // Find n(x) at gaussian quadrature points
      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
        nVec(componentIndex) = nPtr[componentIndex];
      Eigen::VectorXd phiAtQuadPoints = interpMatrix*nVec;
      
      // Compute projection of ln(n_i(x)/nEdge) onto basis functions
      for (int componentIndex = 0; componentIndex < phiAtQuadPoints.rows(); componentIndex++)
      {
        //if (phiAtQuadPoints(componentIndex) < 0.0)
        //  throw Lucee::Except("BoltzmannPhiUpdater::update: n at quad point is negative!");
        phiAtQuadPoints(componentIndex) = gaussWeights[componentIndex]*(meanVtSq*ionMass/elementaryCharge*std::log(phiAtQuadPoints(componentIndex)/nEdge) + phiS);
        //phiAtQuadPoints(componentIndex) = gaussWeights[componentIndex]*(meanVtSq*ionMass/elementaryCharge*std::log(phiAtQuadPoints(componentIndex)/nEdge));
      }

      Eigen::VectorXd phiWeights = massMatrixInv*interpMatrixTranspose*phiAtQuadPoints;

      // Copy into output pointer
      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
        phiPtr[componentIndex] = phiWeights(componentIndex);
    }

    return Lucee::UpdaterStatus();
  }

  void
  BoltzmannPhiUpdater::declareTypes()
  {
    // takes three inputs (moment0=<1>=n(x), vt^2(x), moment1=<v>) 
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // returns one output: mean(vt^2(x)) * ln(n(x)/n(xMax))
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }

  void
  BoltzmannPhiUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}
