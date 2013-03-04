/**
 * @file	LcModalDg1DSymmetricDDGUpdater.cpp
 *
 * @brief	Updater to solver 1D hyperbolic equations using modal DG scheme
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAlignedRectCoordSys.h>
#include <LcMathLib.h>
#include <LcModalDg1DSymmetricDDGUpdater.h>
#include <LcStructuredGridBase.h>

// std includes
#include <vector>

namespace Lucee
{
/** Class id: this is used by registration system */
  const char *ModalDg1DSymmetricDDGUpdater::id = "ModalDg1DSymmetricDDGUpdater";

  ModalDg1DSymmetricDDGUpdater::ModalDg1DSymmetricDDGUpdater()
    : UpdaterIfc()
  {
  }

  void
  ModalDg1DSymmetricDDGUpdater::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::UpdaterIfc::readInput(tbl);
    
    // get number of basis functions to project on
    numBasis = (unsigned) tbl.getNumber("numBasis");

    // Get diffusion coefficient
    diffCoef = (double) tbl.getNumber("diffCoef");

    cfl = tbl.getNumber("cfl"); // CFL number
    cflm = 1.1*cfl; // use slightly large max CFL if not explicitly specified
    if (tbl.hasNumber("cflm"))
      cflm = tbl.getNumber("cflm"); // maximum CFL number

    // check to see if we are only computing L(q) or q+dt*L(q)
    onlyIncrement = false;
    if (tbl.hasBool("onlyIncrement"))
      onlyIncrement = tbl.getBool("onlyIncrement");

    // Check to see if beta0 and beta1 are specified
    if (tbl.hasNumber("beta0"))
    {
      beta0 = tbl.getNumber("beta0");
    }
    else
    {
      if(numBasis == 1)
        beta0 = 0.5;
      else if (numBasis == 2)
        beta0 = 1.5;
      else if (numBasis == 3)
        beta0 = 1.5;
      else
        beta0 = 2;
    }

    if (tbl.hasNumber("beta1"))
    {
      beta1 = tbl.getNumber("beta1");
    }
    else
    {
      if (numBasis < 3)
        beta1 = 0;
      else if (numBasis == 3)
        beta1 = 0.25;
      else
        beta1 = 0.25;
    }
  }

  void
  ModalDg1DSymmetricDDGUpdater::initialize()
  {
    // Get hold of grid
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();
    double dx = grid.getDx(0);
    // Compute weights and ordinates defined on [-1,1] for doing Int[u_x*v_x]
    int integrationDegree = 2*(numBasis-2);
    if (integrationDegree < 0)
      integrationDegree = 0;

    unsigned numGaussPoints = (unsigned)((integrationDegree+1)/2.0 + 0.5);
    std::vector<double> gaussPoints(numGaussPoints);
    std::vector<double> gaussWeights(numGaussPoints);
    legendre_set(numGaussPoints, &gaussPoints[0], &gaussWeights[0]);

    Eigen::MatrixXd DPmk(numBasis, numGaussPoints);

    for (int polyOrder = 0; polyOrder < numBasis; polyOrder++)
    {
      for (int gaussIndex = 0; gaussIndex < numGaussPoints; gaussIndex++)
      {
        // Compute derivatives of Legendre polynomials at ordinates
        DPmk(polyOrder,gaussIndex) = Lucee::legendrePolyDeriv(polyOrder, gaussPoints[gaussIndex]);
      }
    }

    legendreDerivProjections = Eigen::MatrixXd::Zero(numBasis, numBasis);
    double integrationSum;

    // Now evaluate Int[DPj*DPk]
    for (int j = 0; j < numBasis; j ++)
    {
      for (int k = 0; k < numBasis; k++)
      {
        integrationSum = 0.0;
        // Loop over the gauss integration points
        for (int gaussIndex = 0; gaussIndex < numGaussPoints; gaussIndex++)
          integrationSum += gaussWeights[gaussIndex]*DPmk(j,gaussIndex)*DPmk(k,gaussIndex);
        
        legendreDerivProjections(j,k) = (2.0/dx)*integrationSum;
      }
    }

    // Evaluate second derivatives of legendre polynomials at left and right edges
    secondDerivEdgeEvals = Eigen::MatrixXd::Zero(numBasis, 2);

    double sgn = 1.0;

    for (int polyOrder = 0; polyOrder < numBasis; polyOrder++)
    {
      if (polyOrder < 2)
      {
        secondDerivEdgeEvals(polyOrder,0) = 0.0;
        secondDerivEdgeEvals(polyOrder,1) = 0.0;
      }
      else
      {
        secondDerivEdgeEvals(polyOrder,0) = 4.0/(dx*dx)*((2*polyOrder-1)*(sgn*(polyOrder-1)*polyOrder - 
                                              secondDerivEdgeEvals(polyOrder-1,0)) - 
                                              (polyOrder-1)*secondDerivEdgeEvals(polyOrder-2,0))/polyOrder;
        secondDerivEdgeEvals(polyOrder,1) = 4.0/(dx*dx)*((2*polyOrder-1)*((polyOrder-1)*polyOrder + 
                                              secondDerivEdgeEvals(polyOrder-1,1)) - 
                                              (polyOrder-1)*secondDerivEdgeEvals(polyOrder-2,1))/polyOrder;
      }

      sgn *= -1.0;
    }

    // Compute normalization coefficients
    normCoeff = std::vector<double>(numBasis);
    for (unsigned m=0; m<numBasis; ++m)    
      normCoeff[m] = 1/(2.0*m+1);
    
    // Call base class method
    Lucee::UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus
  ModalDg1DSymmetricDDGUpdater::update(double t)
  {
// The algorithms does the first order forward Euler update. It works
// in two stages. In the first stages the increment in the solution,
// i.e. dt*L(q) is computed and stored in qNew. Then, another sweep
// over the domain updates qNew to move the solution in time.

    // Get hold of grid
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();
    // Get input/output arrays
    const Lucee::Field<1, double>& q = this->getInp<Lucee::Field<1, double> >(0);
    Lucee::Field<1, double>& qNew = this->getOut<Lucee::Field<1, double> >(0);

    // Clear out increment
    qNew = 0.0;

    // Time-step
    double dt = t-this->getCurrTime();
    // Get grid size (assume constant size)
    double dx = grid.getDx(0);
    // Figure out whether time step is too large
    if(dt > 1.01*(cfl*dx*dx/(2*diffCoef)))
    {
      // time-step was too large: return a suggestion with correct time-step
      return Lucee::UpdaterStatus(false, cfl*dx*dx/(2*diffCoef));
    }

    // Local region to index
    Lucee::Region<1, int> localRgn = grid.getLocalRegion();

    // state on left and right of each edge
    double qL, qR;
    // first derivative of solution on left and right of each edge
    double qLd1, qRd1;
    // second derivative of solution on left and right of each edge
    double qLd2, qRd2;

    if ( (q.getNumComponents() != numBasis) && (qNew.getNumComponents() != numBasis) )
    {
      Lucee::Except lce(
        "ModalDg1DSymmetricDDGUpdater::update: Number of components in input/output fields should be ");
      lce << numBasis << ". Instead provided " << q.getNumComponents() << std::endl;
    }

    // Iterators
    Lucee::ConstFieldPtr<double> qPtr = q.createConstPtr();
    Lucee::ConstFieldPtr<double> qlPtr = q.createConstPtr();
    Lucee::ConstFieldPtr<double> qrPtr = q.createConstPtr();

    Lucee::FieldPtr<double> qNewPtr = qNew.createPtr();
    Lucee::FieldPtr<double> qNewlPtr = qNew.createPtr();
    Lucee::FieldPtr<double> qNewrPtr = qNew.createPtr();

    // Lower and upper bounds of slice.
    int sliceLower = localRgn.getLower(0);
    int sliceUpper = localRgn.getUpper(0);

    double qFlux;
    double vFluxLeft;
    double vFluxRight;
    int sgn;
// loop over edges computing edge fluxes, accumulating contribution
// from edge flux in cells connected to that edge. NOTE: There is one
// more edge than cells hence the upper limit in the for loop.
    for (int i=sliceLower; i<sliceUpper+1; ++i)
    {
      // Attach iterators to left/right cells of this edge
      q.setPtr(qlPtr, i-1); // left cell
      q.setPtr(qrPtr, i); // right cell

      // Compute conserved variables at left and right of edge
      // Could combine this with below loop...
      evalExpansionRightEdge(qlPtr, qL); // qL is right edge of left cell
      evalExpansionLeftEdge(qrPtr, qR); // qR is left edge of right cell

      // Evaluate the flux in solution at the edge
      // Compute second and first derivative of solution on each side of edge
      qLd1 = 0.0;
      qRd1 = 0.0;
      qLd2 = 0.0;
      qRd2 = 0.0;
      sgn = 1;
      for (int polyOrder = 1; polyOrder < numBasis; polyOrder++)
      {
        qLd1 += polyOrder*(polyOrder+1)/dx*qlPtr[polyOrder];
        qRd1 += sgn*polyOrder*(polyOrder+1)/dx*qrPtr[polyOrder];
        qLd2 += secondDerivEdgeEvals(polyOrder,1)*qlPtr[polyOrder];
        qRd2 += secondDerivEdgeEvals(polyOrder,0)*qrPtr[polyOrder];
        sgn *= -1;
      }
      qFlux = beta0*(qR-qL)/dx + 0.5*(qLd1+qRd1) + beta1*dx*(qRd2-qLd2);

      // Attach iterators to left and right cells to accumulate increment
      qNew.setPtr(qNewlPtr, i-1); // left cell
      qNew.setPtr(qNewrPtr, i); // right cell

      // accumulate increment to appropriate cells
      sgn = 1;
      for (int m = 0; m < numBasis; m++)
      {
        // vFluxLeft corresponds to interface term _{j-1/2}, but named because this
        // contribution will be added to the left cell
        vFluxLeft = (qR-qL)*(-beta0/dx + m*(m+1)/(2.0*dx) - beta1*dx*secondDerivEdgeEvals(m,1));
        vFluxRight = (qR-qL)*(sgn*beta0/dx - sgn*m*(m+1)/(2.0*dx) + beta1*dx*secondDerivEdgeEvals(m,0));
        qNewlPtr[m] += qFlux - vFluxLeft;
        qNewrPtr[m] += -sgn*qFlux - vFluxRight;
        sgn *= -1;
      }
    }
    
    double integralResult = 0.0;
    // loop over cells adding contribution from volume integrals
    // This will only be needed if numBasis > 1
    if (numBasis > 1)
    {
      for (unsigned i=sliceLower; i<sliceUpper; ++i)
      {
        // Attach iterators to cell
        q.setPtr(qPtr, i);
        qNew.setPtr(qNewPtr, i);
        // Loop over each basis function. Only when test function is
        // degree 1 or higher will contribute due to (v_m)_{x} term
        for (int basisIndex = 1; basisIndex < numBasis; basisIndex++)
        {
          integralResult = 0.0;
          // Compute Integral[(v_m)_{x}*u_{x} dx] over cell by summing
          // weights of u with precomputed integral
          for (int uCoeffIndex = 0; uCoeffIndex < numBasis; uCoeffIndex++)
          {
            integralResult += qPtr[uCoeffIndex]*legendreDerivProjections(basisIndex,uCoeffIndex);
          }
          // Add this result to the new state of the cell
          qNewPtr[basisIndex] -= integralResult;
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
  ModalDg1DSymmetricDDGUpdater::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }

  void
  ModalDg1DSymmetricDDGUpdater::evalExpansionLeftEdge(const Lucee::ConstFieldPtr<double>& qCoeff,
    double& qOut)
  {
    int sgn = 1.0;
    qOut = 0.0;
    for (unsigned m=0; m<numBasis; ++m)
    {
      qOut += sgn*qCoeff[m];
      sgn *= -1;
    }
  }

  void
  ModalDg1DSymmetricDDGUpdater::evalExpansionRightEdge(const Lucee::ConstFieldPtr<double>& qCoeff,
    double& qOut)
  {
    qOut = 0.0;
    for (unsigned m=0; m<numBasis; ++m)
    {
      qOut += qCoeff[m];
    }
  }
}
