/**
 * @file	LcModalDg1DLocalDGUpdater.cpp
 *
 * @brief	Updater to solver 1D diffusion equations using the local dg scheme.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAlignedRectCoordSys.h>
#include <LcMathLib.h>
#include <LcModalDg1DLocalDGUpdater.h>
#include <LcStructuredGridBase.h>

// std includes
#include <vector>

namespace Lucee
{
/** Class id: this is used by registration system */
  const char *ModalDg1DLocalDGUpdater::id = "ModalDg1DLocalDG";

  ModalDg1DLocalDGUpdater::ModalDg1DLocalDGUpdater()
    : UpdaterIfc()
  {
  }

  void
  ModalDg1DLocalDGUpdater::readInput(Lucee::LuaTable& tbl)
  {
    // Call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    // Get number of basis functions to project on
    numBasis = (unsigned) tbl.getNumber("numBasis");

    // Get diffusion coefficient
    diffCoef = (double) tbl.getNumber("diffCoef");
    // Actually sqrt it
    diffCoefSqrt = std::sqrt(diffCoef);

    cfl = tbl.getNumber("cfl"); // CFL number
    cflm = 1.1*cfl; // use slightly large max CFL if not explicitly specified
    if (tbl.hasNumber("cflm"))
      cflm = tbl.getNumber("cflm"); // maximum CFL number

    fluxPair = AVERAGE;
    if (tbl.hasString("fluxPair"))
    {
      if (tbl.getString("fluxPair") == "left")
        fluxPair = LEFT;
      else if (tbl.getString("fluxPair") == "right")
        fluxPair = RIGHT;
      else if (tbl.getString("fluxPair") == "average")
        fluxPair = AVERAGE;
      else
      {
        Lucee::Except lce("LcModalDg1DLocalDGUpdater::readInput: 'fluxPair' must be one of ");
        lce << " 'left' 'right' 'average'. Provided '" << tbl.getString("fluxPair")
            << "' instead";
        throw lce;
      }
    }

    onlyIncrement = false;
    if (tbl.hasBool("onlyIncrement"))
      onlyIncrement = tbl.getBool("onlyIncrement");
  }


  void
  ModalDg1DLocalDGUpdater::initialize()
  {
    // Get hold of grid
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();
    double dx = grid.getDx(0);

    // Precompute volume integrals if they will be needed
    if(numBasis > 1)
    {
      // For computation of P_m*P_n' integrals
      int integrationDegree = 2*numBasis - 3;
      // Compute weights and ordinates defined on [-1,1] for doing recovery integrals
      unsigned numGaussPoints = (unsigned)((integrationDegree+1)/2.0 + 0.5);
      std::vector<double> gaussPoints(numGaussPoints);
      std::vector<double> gaussWeights(numGaussPoints);
      legendre_set(numGaussPoints, &gaussPoints[0], &gaussWeights[0]);
      // Compute elements of P_m' P_n matrix
      volumeMatrix = Eigen::MatrixXd::Zero(numBasis,numBasis);
      double integralValue;

      for(int rowIndex = 0; rowIndex < volumeMatrix.rows(); rowIndex++)
      {
        // Fill out rows 2*rowIndex and 2*rowIndex+1
        for(int colIndex = 0; colIndex < volumeMatrix.cols(); colIndex++)
        {
          integralValue = 0.0;
          // Compute integral of P_m' * P_n using gaussian quadrature
          for(int nodeIndex = 0; nodeIndex < numGaussPoints; nodeIndex++)
          {
            integralValue += gaussWeights[nodeIndex]*
              Lucee::legendrePolyDeriv(rowIndex,gaussPoints[nodeIndex])*
              Lucee::legendrePoly(colIndex,gaussPoints[nodeIndex]);
          }
          // Store results in the matrix
          volumeMatrix(rowIndex,colIndex) = integralValue;
        }
      }
    }

    normCoeff = std::vector<double>(numBasis);
    // Compute normalization coefficients
    for (unsigned m=0; m<numBasis; ++m)    
      normCoeff[m] = -1.0*diffCoefSqrt*(2.0*m+1)/dx;
    // Call base class method
    Lucee::UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus
  ModalDg1DLocalDGUpdater::update(double t)
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

    // wPlus represents using flux f = f+
    Lucee::Field<1, double> wPlus = qNew.duplicate();
    // wMinus represents using flux f = f-
    Lucee::Field<1, double> wMinus = qNew.duplicate();
    // Clear out w
    wPlus = 0.0;
    wMinus = 0.0;

    // clear out increment
    qNew = 0.0;

    // time-step
    double dt = t-this->getCurrTime();
    double dx = grid.getDx(0);

    // Figure out whether time step is too large
    // TODO: determine whether or not to keep 1.01 factor
    if(dt > 1.01*(cfl*dx*dx/(2*diffCoef)))
    {
      // time-step was too large: return a suggestion with correct time-step
      return Lucee::UpdaterStatus(false, cfl*dx*dx/(2*diffCoef));
    }

    // local region to index
    Lucee::Region<1, int> localRgn = grid.getLocalRegion();

    if ( (q.getNumComponents() != numBasis) && (qNew.getNumComponents() != numBasis) )
    {
      Lucee::Except lce(
        "ModalDg1DLocalDGUpdater::update: Number of components in input/output fields should be ");
      lce << numBasis << ". Instead provided " << q.getNumComponents() << std::endl;
    }

    // iterators
    Lucee::ConstFieldPtr<double> qPtr = q.createConstPtr();
    Lucee::ConstFieldPtr<double> qlPtr = q.createConstPtr();
    Lucee::ConstFieldPtr<double> qrPtr = q.createConstPtr();

    Lucee::FieldPtr<double> qNewPtr = qNew.createPtr();
    Lucee::FieldPtr<double> qNewlPtr = qNew.createPtr();
    Lucee::FieldPtr<double> qNewrPtr = qNew.createPtr();
    
    Lucee::FieldPtr<double> wPlusPtr = wPlus.createPtr();
    Lucee::FieldPtr<double> wPluslPtr = wPlus.createPtr();
    Lucee::FieldPtr<double> wPlusrPtr = wPlus.createPtr();
    
    Lucee::FieldPtr<double> wMinusPtr = wMinus.createPtr();
    Lucee::FieldPtr<double> wMinuslPtr = wMinus.createPtr();
    Lucee::FieldPtr<double> wMinusrPtr = wMinus.createPtr();

    // lower and upper bounds of slice.
    int sliceLower = localRgn.getLower(0);
    int sliceUpper = localRgn.getUpper(0);

    // loop over edges computing edge fluxes, accumulating contribution
    // from edge flux in cells connected to that edge. NOTE: There is one
    // more edge than cells hence the upper limit in the for loop.
    double wRight;
    double wLeft;
    int sgn = -1;

    // First compute the updated w_i's
    for (int i=sliceLower; i<sliceUpper+1; ++i)
    {
      // Attach iterators to left/right cells of this edge
      q.setPtr(qlPtr, i-1); // left cell
      q.setPtr(qrPtr, i); // right cell
      wPlus.setPtr(wPluslPtr, i-1);
      wPlus.setPtr(wPlusrPtr, i);
      wMinus.setPtr(wMinuslPtr, i-1);
      wMinus.setPtr(wMinusrPtr, i);
      
      wRight = 0.0;
      wLeft = 0.0;
      sgn = 1;

      for (int m = 0; m < numBasis; m++)
      {
        wRight += sgn*qrPtr[m];
        wLeft  += qlPtr[m];
        sgn *= -1;
      }

      // Accumulate increment to appropriate cells
      sgn = -1;
      for (int m = 0; m < numBasis; m++)
      {
        wPluslPtr[m] += wRight;
        wPlusrPtr[m] += sgn*wRight;
        wMinuslPtr[m] += wLeft;
        wMinusrPtr[m] += sgn*wLeft;
        sgn *= -1;
      }
    }

    double integralResult;
    // Loop over cells adding contribution from volume integrals
    // This will only be needed if numBasis > 1
    if (numBasis > 1)
    {
      for (int i = sliceLower; i < sliceUpper; i++)
      {
        // Attach iterators to cell
        q.setPtr(qPtr, i);
        wPlus.setPtr(wPlusPtr, i);
        wMinus.setPtr(wMinusPtr, i);
        // Loop over each basis function. Only when test function is
        // degree 1 or higher will contribute due to (v_m)_{x} term
        for (int basisIndex = 1; basisIndex < numBasis; basisIndex++)
        {
          integralResult = 0.0;
          // Compute Integral[(v_m)_{x}*w dx] over cell by summing
          // weights of w with precomputed integral
          for (int wCoeffIndex = 0; wCoeffIndex < numBasis; wCoeffIndex++)
          {
            integralResult += qPtr[wCoeffIndex]*volumeMatrix(basisIndex,wCoeffIndex);
          }
          // Add this result to the new state of the cell
          wPlusPtr[basisIndex] -= integralResult;
          wMinusPtr[basisIndex] -= integralResult;
        }
      }
    }
    
    for (int i = sliceLower; i < sliceUpper; i++)
    {
      wPlus.setPtr(wPlusPtr, i);
      wMinus.setPtr(wMinusPtr, i);
      // Normalize increment
      for (int m = 0; m < numBasis; m++)
      {
        wPlusPtr[m] *= normCoeff[m];
        wMinusPtr[m] *= normCoeff[m];
      }
    }

    wPlus.applyPeriodicBc(0);
    wMinus.applyPeriodicBc(0);
    
    double qLeft;
    double qRight;

    // Now use updatedWRight to compute g = df/dt update
    for (int i=sliceLower; i<sliceUpper+1; ++i)
    {
      // Attach iterators to left/right cells of this edge
      wPlus.setPtr(wPluslPtr, i-1);
      wPlus.setPtr(wPlusrPtr, i);
      wMinus.setPtr(wMinuslPtr, i-1);
      wMinus.setPtr(wMinusrPtr, i);
      // Attach iterators to left and right cells to accumulate increment
      qNew.setPtr(qNewlPtr, i-1); // left cell
      qNew.setPtr(qNewrPtr, i); // right cell
      
      qLeft = 0.0;
      qRight = 0.0;
      sgn = 1;

      for (int m = 0; m < numBasis; m++)
      {
        qLeft  += wPluslPtr[m];
        qRight += sgn*wMinusrPtr[m];
        sgn *= -1;
      }

      // Accumulate increment to appropriate cells
      sgn = -1;
      for (int m = 0; m < numBasis; m++)
      {
        if (fluxPair == RIGHT)
        {
          qNewlPtr[m] += qLeft;
          qNewrPtr[m] += sgn*qLeft;
        }
        else if (fluxPair == LEFT)
        {
          qNewlPtr[m] += qRight;
          qNewrPtr[m] += sgn*qRight;
        }
        else
        {
          qNewlPtr[m] += 0.5*(qRight+qLeft);
          qNewrPtr[m] += sgn*0.5*(qRight+qLeft);
        }
        
        sgn *= -1;
      }
    }

    // Loop over cells adding contribution from volume integrals
    // This will only be needed if numBasis > 1
    if (numBasis > 1)
    {
      for (int i = sliceLower; i < sliceUpper; i++)
      {
        // Attach iterators to cell
        qNew.setPtr(qNewPtr, i);
        wPlus.setPtr(wPlusPtr, i);
        wMinus.setPtr(wMinusPtr, i);
        // Loop over each basis function. Only when test function is
        // degree 1 or higher will contribute due to (v_m)_{x} term
        for (int basisIndex = 1; basisIndex < numBasis; basisIndex++)
        {
          integralResult = 0.0;
          // Compute Integral[(v_m)_{x}*w dx] over cell by summing
          // weights of w with precomputed integral
          for (int wCoeffIndex = 0; wCoeffIndex < numBasis; wCoeffIndex++)
          {
            if (fluxPair == RIGHT)
            {
              integralResult += wPlusPtr[wCoeffIndex]*
                volumeMatrix(basisIndex,wCoeffIndex);
            }
            else if (fluxPair == LEFT)
            {
              integralResult += wMinusPtr[wCoeffIndex]*
                volumeMatrix(basisIndex,wCoeffIndex);
            }
            else
            {
              integralResult += 0.5*(wMinusPtr[wCoeffIndex] + wPlusPtr[wCoeffIndex])*
                volumeMatrix(basisIndex,wCoeffIndex);
            }
            
          }
          // Add this result to the new state of the cell
          qNewPtr[basisIndex] -= integralResult;
        }
      }
    }

    // Perform final update
    for (int i = sliceLower; i < sliceUpper; i++)
    {
      qNew.setPtr(qNewPtr, i);
      // Normalize increment
      for (int m = 0; m < numBasis; m++)
      {
        qNewPtr[m] *= normCoeff[m];
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
  ModalDg1DLocalDGUpdater::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }
}
