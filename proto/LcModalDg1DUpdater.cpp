/**
 * @file	LcModalDg1DUpdater.cpp
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
#include <LcModalDg1DUpdater.h>
#include <LcStructuredGridBase.h>

namespace Lucee
{
/** Class id: this is used by registration system */
  const char *ModalDg1DUpdater::id = "ModalDg1D";

  ModalDg1DUpdater::ModalDg1DUpdater()
    : UpdaterIfc(), numBasis(1), Pmk(1,1), DPmk(1,1), normCoeff(1), w(1), mu(1)
  {
  }

  void
  ModalDg1DUpdater::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::UpdaterIfc::readInput(tbl);

// equation to solve
    if (tbl.hasObject<Lucee::HyperEquation>("equation"))
      equation = &tbl.getObjectAsBase<Lucee::HyperEquation>("equation");
    else
    {
      Lucee::Except lce("ModalDg1DUpdater::readInput: Must specify an equation to solve!");
      throw lce;
    }
    meqn = equation->getNumEqns();

// get number of basis functions to project on
    numBasis = (unsigned) tbl.getNumber("numBasis");

    cfl = tbl.getNumber("cfl"); // CFL number
    cflm = 1.1*cfl; // use slightly large max CFL if not explicitly specified
    if (tbl.hasNumber("cflm"))
      cflm = tbl.getNumber("cflm"); // maximum CFL number

// allocate space
    Pmk = Matrix<double>(numBasis, numBasis);
    DPmk = Matrix<double>(numBasis, numBasis);
    normCoeff = Lucee::Vector<double>(numBasis);
    w = Lucee::Vector<double>(numBasis);
    mu = Lucee::Vector<double>(numBasis);

// compute weights and ordinates
    Lucee::gauleg(numBasis, -1, 1, mu, w);

// compute Legendre polynomials at ordinates
    for (unsigned m=0; m<numBasis; ++m)
      for (unsigned k=0; k<numBasis; ++k)
        Pmk(m,k) = Lucee::legendrePoly(m, mu[k]);

// compute derivatives of Legendre polynomials at ordinates
    for (unsigned m=0; m<numBasis; ++m)
      for (unsigned k=0; k<numBasis; ++k)
        DPmk(m,k) = Lucee::legendrePolyDeriv(m, mu[k]);

// compute normalization coefficients
    for (unsigned m=0; m<numBasis; ++m)    
      normCoeff[m] = 1/(2.0*m+1);
  }

  void
  ModalDg1DUpdater::initialize()
  {
// call base class method
    Lucee::UpdaterIfc::initialize();

  }

  Lucee::UpdaterStatus
  ModalDg1DUpdater::update(double t)
  {
// get hold of grid
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();
// get input/output arrays
    const Lucee::Field<1, double>& q = this->getInp<Lucee::Field<1, double> >(0);
    Lucee::Field<1, double>& qNew = this->getOut<Lucee::Field<1, double> >(0);

// clear out output field
    qNew = 0.0;

// time-step
    double dt = t-this->getCurrTime();
// local region to index
    Lucee::Region<1, int> localRgn = grid.getLocalRegion();

// create coordinate system along X-axis
    Lucee::AlignedRectCoordSys coordSys(0);

// state on left and right of each edge
    Lucee::FieldPtr<double> qL(meqn), qR(meqn);
// flux and numerical flux
    Lucee::FieldPtr<double> flux(meqn), numFlux(meqn);

    if ( (q.getNumComponents() != meqn*numBasis) && (qNew.getNumComponents() != meqn*numBasis) )
    {
      Lucee::Except lce(
        "ModalDg1DUpdater::update: Number of components in input/output fields should be ");
      lce << meqn*numBasis << ". Instead provided " << q.getNumComponents() << std::endl;
    }

// iterators to cells on left/right of edge
    Lucee::ConstFieldPtr<double> qPtr = q.createConstPtr();
    Lucee::ConstFieldPtr<double> qlPtr = q.createConstPtr();
    Lucee::ConstFieldPtr<double> qrPtr = q.createConstPtr();

    Lucee::FieldPtr<double> qNewPtr = qNew.createPtr();
    Lucee::FieldPtr<double> qNewlPtr = qNew.createPtr();
    Lucee::FieldPtr<double> qNewrPtr = qNew.createPtr();

// lower and upper bounds of slice.
    int sliceLower = localRgn.getLower(0);
    int sliceUpper = localRgn.getUpper(0);

    double dxL, dxR, dx, xc[3];
// maximum CFL number used
    double cfla = 0.0;

// The algorithms does the first order forward Euler update. It works
// in two stages. In the first stages the increment in the solution,
// i.e. dt*L(q) is computed and stored in qNew. Then, another sweep
// over the domain updates qNew = q + dt*L(q) to move the solution in
// time.

// loop over edges computing edge fluxes, accumulating contribution
// from edge flux in cells connected to that edge. NOTE: There is one
// more edge that cells hence the upper limit in the for loop.
    for (int i=sliceLower; i<sliceUpper+1; ++i)
    {
// cell spacing in left cell
      grid.setIndex(i-1);
      dxL = grid.getDx(0);

// cell spacing in right cell
      grid.setIndex(i);
      dxR = grid.getDx(0);

// attach iterators to left/right cells of this edge
      q.setPtr(qlPtr, i-1); // left cell
      q.setPtr(qrPtr, i); // right cell

// compute conserved variables at left and right of edge
      evalExpansionRightEdge(qlPtr, qL); // qL is right edge of left cell
      evalExpansionLeftEdge(qrPtr, qR); // qR is left edge of right cell

// compute numerical flux at edge
      double maxs = equation->numericalFlux(coordSys, qL, qR, numFlux);

// attach iterators to left and right cells to accumulate increment
      qNew.setPtr(qNewlPtr, i-1); // left cell
      qNew.setPtr(qNewrPtr, i); // right cell

// accumulate increment to appropriate cells
      for (unsigned k=0; k<meqn; ++k)
      {
        int sgn = 1.0;
        for (unsigned m=0; m<numBasis; ++m)
        {
          qNewlPtr[k+m*meqn] -= numFlux[k];
          qNewrPtr[k+m*meqn] += sgn*numFlux[k];
          sgn *= -1;
        }
      }

      double dtdx = 2*dt/(dxL+dxR);
// compute current maximum CFL number
      cfla = Lucee::max3(cfla, dtdx*maxs, -dtdx*maxs);
    }

    if (cfla > cflm)
// time-step was too large: return a suggestion with correct time-step
      return Lucee::UpdaterStatus(false, dt*cfl/cfla);

// loop over cells adding contribution from volume integrals
    for (unsigned i=sliceLower; i<sliceUpper; ++i)
    {
// attach iterators to cell
      q.setPtr(qPtr, i);
      qNew.setPtr(qNewPtr, i);

// loop over each ordinate computing flux
      for (unsigned r=0; r<numBasis; ++r)
      {
// compute conserved variable at ordinate and store in qL
        evalExpansion(qPtr, r, qL);
// compute flux at this location
        equation->flux(coordSys, qL, flux);

// accumulate contribution of this flux to each mode for each equation
        for (unsigned m=1; m<numBasis; ++m)
        { // no contribution to m=0 mode
          for (unsigned k=0; k<meqn; ++k)
          {
            qNewPtr[k+m*meqn] += w[r]*flux[k]*DPmk(m,r);
          }
        }
      }
    }

// perform final update
    double nc;
    for (unsigned i=sliceLower; i<sliceUpper; ++i)
    {
      grid.setIndex(i);
      double dx = grid.getDx(0);

      qNew.setPtr(qNewPtr, i);
// normalize increment
      for (unsigned m=0; m<numBasis; ++m)
      {
        nc = dt/(normCoeff[m]*dx);
        for (unsigned k=0; k<meqn; ++k)
          qNewPtr[k+m*meqn] *= nc;
      }

      q.setPtr(qPtr, i);
// do forward Euler. NOTE: qNew at this point has dt*L(q) and after
// will have q + dt*L(q).
      for (unsigned n=0; n<q.getNumComponents(); ++n)
        qNewPtr[n] += qPtr[n];
    }

    return Lucee::UpdaterStatus(true, dt*cfl/cfla);
  }

  void
  ModalDg1DUpdater::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }

  void
  ModalDg1DUpdater::evalExpansion(const Lucee::ConstFieldPtr<double>& qCoeff,
    unsigned c, Lucee::FieldPtr<double>& qOut)
  {
    for (unsigned k=0; k<meqn; ++k)
    {
      qOut[k] = 0.0;
      for (unsigned m=0; m<numBasis; ++m)
        qOut[k] += qCoeff[k+m*meqn]*Pmk(m, c);
    }
  }

  void
  ModalDg1DUpdater::evalExpansionLeftEdge(const Lucee::ConstFieldPtr<double>& qCoeff,
    Lucee::FieldPtr<double>& qOut)
  {
    for (unsigned k=0; k<meqn; ++k)
    {
      int sgn = 1.0;
      qOut[k] = 0.0;
      for (unsigned m=0; m<numBasis; ++m)
      {
        qOut[k] += sgn*qCoeff[k+m*meqn];
        sgn *= -1;
      }
    }
  }

  void
  ModalDg1DUpdater::evalExpansionRightEdge(const Lucee::ConstFieldPtr<double>& qCoeff,
    Lucee::FieldPtr<double>& qOut)
  {
    for (unsigned k=0; k<meqn; ++k)
    {
      qOut[k] = 0.0;
      for (unsigned m=0; m<numBasis; ++m)
      {
        qOut[k] += qCoeff[k+m*meqn];
      }
    }
  }
}
