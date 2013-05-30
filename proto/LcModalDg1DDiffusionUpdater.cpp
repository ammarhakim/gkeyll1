/**
 * @file	LcModalDg1DDiffusionUpdater.cpp
 *
 * @brief	Updater to solver 1D diffusion equations using modal DG scheme
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAlignedRectCoordSys.h>
#include <LcMathLib.h>
#include <LcModalDg1DDiffusionUpdater.h>
#include <LcStructuredGridBase.h>

// std includes
#include <vector>

namespace Lucee
{
/** Class id: this is used by registration system */
  const char *ModalDg1DDiffusionUpdater::id = "ModalDg1DDiffusion";

  ModalDg1DDiffusionUpdater::ModalDg1DDiffusionUpdater()
    : UpdaterIfc()
  {
  }

  void
  ModalDg1DDiffusionUpdater::readInput(Lucee::LuaTable& tbl)
  {
    // Call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    // Get number of basis functions to project on
    numBasis = (unsigned) tbl.getNumber("numBasis");

    // Get diffusion coefficient
    diffCoef = (double) tbl.getNumber("diffCoef");

    cfl = tbl.getNumber("cfl"); // CFL number
    cflm = 1.1*cfl; // use slightly large max CFL if not explicitly specified
    if (tbl.hasNumber("cflm"))
      cflm = tbl.getNumber("cflm"); // maximum CFL number

    onlyIncrement = false;
    if (tbl.hasBool("onlyIncrement"))
      onlyIncrement = tbl.getBool("onlyIncrement");
  }


  void
  ModalDg1DDiffusionUpdater::initialize()
  {
    // Get hold of grid
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();
    double dx = grid.getDx(0);
    // Assume polyOrder = numBasis - 1
    int recoveryDegree = (2*numBasis-1);
    // For computation of legendreP*x^n integarls
    int integrationDegree = (numBasis-1)+recoveryDegree;

    // Representation of recovery polynomial
    std::vector<double> recoveryPoly(recoveryDegree+1);

    // Fill out coefficients of recovery polynomial terms
    recoveryPoly[0] = 1.0;
    for(int termIndex = 1; termIndex < recoveryPoly.size(); termIndex++)
    {
      recoveryPoly[termIndex] = recoveryPoly[termIndex-1]/termIndex;
    }

    // Compute weights and ordinates defined on [-1,1] for doing recovery integrals
    unsigned numGaussPoints = (unsigned)((integrationDegree+1)/2.0 + 0.5);
    std::vector<double> gaussPoints(numGaussPoints);
    std::vector<double> gaussWeights(numGaussPoints);
    legendre_set(numGaussPoints, &gaussPoints[0], &gaussWeights[0]);

    Eigen::MatrixXd legendreEvaluations(numGaussPoints,numBasis);
    // Evaluate legendre polynomials at gaussian quadrature points on reference cell
    for(int basisIndex = 0; basisIndex < numBasis; basisIndex++)
    {
      for(int nodeIndex = 0; nodeIndex < numGaussPoints; nodeIndex++)
      {
        legendreEvaluations(nodeIndex,basisIndex) = Lucee::legendrePoly(basisIndex,gaussPoints[nodeIndex]);
      }
    }

    // Compute elements of recovery matrix
    Eigen::MatrixXd recoveryMatrix(recoveryDegree+1,recoveryDegree+1);
    double integralValueLeft;
    double integralValueRight;

    for(int rowIndex = 0; rowIndex < recoveryMatrix.rows()/2; rowIndex++)
    {
      // Fill out rows 2*rowIndex and 2*rowIndex+1
      for(int colIndex = 0; colIndex < recoveryMatrix.cols(); colIndex++)
      {
        integralValueLeft = 0.0;
        integralValueRight = 0.0;
        // Compute integral of legendre polynomial with a power of x using gaussian quadrature
        for(int nodeIndex = 0; nodeIndex < numGaussPoints; nodeIndex++)
        {
          // left cell integral
          integralValueLeft += gaussWeights[nodeIndex]*recoveryPoly[colIndex]*pow((dx/2.0)*(gaussPoints[nodeIndex]-1),colIndex)*
            legendreEvaluations(nodeIndex,rowIndex);

          // right cell integral
          integralValueRight += gaussWeights[nodeIndex]*recoveryPoly[colIndex]*pow((dx/2.0)*(gaussPoints[nodeIndex]+1),colIndex)*
            legendreEvaluations(nodeIndex,rowIndex);
        }
        // Store results in recovery matrix
        recoveryMatrix(2*rowIndex,colIndex) = ((2*rowIndex+1)/2.0)*integralValueLeft;
        recoveryMatrix(2*rowIndex+1,colIndex) = ((2*rowIndex+1)/2.0)*integralValueRight;
      }
    }

    recoveryMatrixInv = recoveryMatrix.inverse();

    // Precompute volume integrals if they will be needed
    if(numBasis > 2)
    {
      double integralValue;
      // Recompute number of gaussian integration points needed
      integrationDegree = 2*numBasis - 4;
      numGaussPoints    = (unsigned)((integrationDegree+1)/2.0 + 0.5);
      gaussPoints.resize(numGaussPoints);
      gaussWeights.resize(numGaussPoints);
      legendre_set(numGaussPoints, &gaussPoints[0], &gaussWeights[0]);

      // Allocate storage matrix
      recoveryVolume = Eigen::MatrixXd::Zero(numBasis,numBasis);
      
      // Start at 2 since we know that first two rows will be zeros
      for(int rowIndex = 2; rowIndex < numBasis; rowIndex++)
      {
        for(int colIndex = 0; colIndex < numBasis; colIndex++)
        {
          integralValue = 0.0;
          // Evalute integral
          for(int nodeIndex = 0; nodeIndex < numGaussPoints; nodeIndex++)
          {
            integralValue += gaussWeights[nodeIndex]*legendrePoly2ndDeriv(rowIndex,gaussPoints[nodeIndex])*
              Lucee::legendrePoly(colIndex,gaussPoints[nodeIndex]);
          }
          recoveryVolume(rowIndex,colIndex) = (2.0/dx)*integralValue;
        }
      }
    }

    normCoeff = std::vector<double>(numBasis);
    // Compute normalization coefficients
    for (unsigned m=0; m<numBasis; ++m)    
      normCoeff[m] = 1/(2.0*m+1);
// call base class method
    Lucee::UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus
  ModalDg1DDiffusionUpdater::update(double t)
  {
// The algorithms does the first order forward Euler update. It works
// in two stages. In the first stages the increment in the solution,
// i.e. dt*L(q) is computed and stored in qNew. Then, another sweep
// over the domain updates qNew to move the solution in time.

// get hold of grid
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();
// get input/output arrays
    const Lucee::Field<1, double>& q = this->getInp<Lucee::Field<1, double> >(0);
    Lucee::Field<1, double>& qNew = this->getOut<Lucee::Field<1, double> >(0);

// clear out increment
    qNew = 0.0;

// time-step
    double dt = t-this->getCurrTime();
    double dx = grid.getDx(0);

    // Figure out whether time step is too large
    if(dt > (1.01*cfl*dx*dx/(2*diffCoef)))
    {
      // time-step was too large: return a suggestion with correct time-step
      return Lucee::UpdaterStatus(false, cfl*dx*dx/(2*diffCoef));
    }

// local region to index
    Lucee::Region<1, int> localRgn = grid.getLocalRegion();

// state on left and right of each edge
    Lucee::FieldPtr<double> qL(1), qR(1);

    if ( (q.getNumComponents() != numBasis) && (qNew.getNumComponents() != numBasis) )
    {
      Lucee::Except lce(
        "ModalDg1DDiffusionUpdater::update: Number of components in input/output fields should be ");
      lce << numBasis << ". Instead provided " << q.getNumComponents() << std::endl;
    }

// iterators
    Lucee::ConstFieldPtr<double> qPtr = q.createConstPtr();
    Lucee::ConstFieldPtr<double> qlPtr = q.createConstPtr();
    Lucee::ConstFieldPtr<double> qrPtr = q.createConstPtr();

    Lucee::FieldPtr<double> qNewPtr = qNew.createPtr();
    Lucee::FieldPtr<double> qNewlPtr = qNew.createPtr();
    Lucee::FieldPtr<double> qNewrPtr = qNew.createPtr();

// lower and upper bounds of slice.
    int sliceLower = localRgn.getLower(0);
    int sliceUpper = localRgn.getUpper(0);

    std::vector<double> boundaryTerms(2);

// loop over edges computing edge fluxes, accumulating contribution
// from edge flux in cells connected to that edge. NOTE: There is one
// more edge that cells hence the upper limit in the for loop.
    for (int i=sliceLower; i<sliceUpper+1; ++i)
    {
// attach iterators to left/right cells of this edge
      q.setPtr(qlPtr, i-1); // left cell
      q.setPtr(qrPtr, i); // right cell

      // Calculate recovery polynomial coefficients f0 and f1
      computeRecoveryPolynomial(qlPtr, qrPtr, boundaryTerms);

// attach iterators to left and right cells to accumulate increment
      qNew.setPtr(qNewlPtr, i-1); // left cell
      qNew.setPtr(qNewrPtr, i); // right cell

// accumulate increment to appropriate cells
      int sgn = -1.0;
      for (int m = 0; m < numBasis; m++)
      {
        qNewlPtr[m] += boundaryTerms[1] - (m*(m+1)/dx)*boundaryTerms[0];
        qNewrPtr[m] += sgn*(boundaryTerms[1] + (m*(m+1)/dx)*boundaryTerms[0]);
        sgn *= -1;
      }
    }

    double integralResult = 0.0;
// loop over cells adding contribution from volume integrals
// This will only be needed if numBasis > 2
    if (numBasis > 2)
    {
      for (unsigned i=sliceLower; i<sliceUpper; ++i)
      {
        // Attach iterators to cell
        q.setPtr(qPtr, i);
        qNew.setPtr(qNewPtr, i);
        // Loop over each basis function. Only when test function is
        // degree 2 or higher will contribute due to (v_m)_{xx} term
        for (int basisIndex = 2; basisIndex < numBasis; basisIndex++)
        {
          integralResult = 0.0;
          // Compute Integral[(v_m)_{xx}*u dx] over cell by summing
          // weights of u with precomputed integral
          for (int uCoeffIndex = 0; uCoeffIndex < numBasis; uCoeffIndex++)
          {
            integralResult += qPtr[uCoeffIndex]*recoveryVolume(basisIndex,uCoeffIndex);
          }
          // Add this result to the new state of the cell
          qNewPtr[basisIndex] += integralResult;
        }
      }
    }

    // Perform final update
    double nc;
    for (unsigned i=sliceLower; i<sliceUpper; ++i)
    {
      qNew.setPtr(qNewPtr, i);
      // Normalize increment
      for (int m = 0; m < numBasis; m++)
      {
        nc = diffCoef/(normCoeff[m]*dx);
        qNewPtr[m] *= nc;
      }

      // Only do this if full update is requested
      if (onlyIncrement == false)
      {
        q.setPtr(qPtr, i);
        // do forward Euler step
        for (int n = 0; n < q.getNumComponents(); n++)
          qNewPtr[n] = qPtr[n] + dt*qNewPtr[n];
      }
    }

    return Lucee::UpdaterStatus(true, cfl*dx*dx/(2*diffCoef));
  }

  void
  ModalDg1DDiffusionUpdater::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }

  double ModalDg1DDiffusionUpdater::legendrePoly2ndDeriv(int order, double x)
  {
    return (-2*x/(x*x-1)*Lucee::legendrePolyDeriv(order, x) + 
      order/(x*x-1)*(x*Lucee::legendrePolyDeriv(order, x) + Lucee::legendrePoly(order, x) -
      Lucee::legendrePolyDeriv(order-1,x)));
  }

  void
  ModalDg1DDiffusionUpdater::computeRecoveryPolynomial(const Lucee::ConstFieldPtr<double>& qLeft,
    Lucee::ConstFieldPtr<double>& qRight, std::vector<double>& boundaryTerms)
  {

    Eigen::MatrixXd qVals(2*numBasis,1);

    // Fill out qVals
    for (int componentIndex = 0; componentIndex < numBasis; componentIndex++)
    {
      // left goes in 2*componentIndex, right goes in 2*componentIndex+1
      qVals(2*componentIndex,0)   = qLeft[componentIndex];
      qVals(2*componentIndex+1,0) = qRight[componentIndex];
    }

    // Get coefficients of reconstructed polynomial at a single edge
    Eigen::MatrixXd fVals(2*numBasis,1);

    fVals = recoveryMatrixInv*qVals;

    // Copy relevant data to boundaryTerms matrix
    boundaryTerms[0] = fVals(0,0);
    boundaryTerms[1] = fVals(1,0);
  }
}
