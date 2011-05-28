/**
 * @file	LcMusclHancock1DUpdater.h
 *
 * @brief	Solver for 1D Euler equations using MUSCl-Hancock scheme.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcMathLib.h>
#include <LcMusclHancock1DUpdater.h>
#include <LcStructuredGridBase.h>

// std includes
#include <cmath>

namespace Lucee
{
// limiter tags
  static const unsigned AVERAGE_LIMITER = 0;
  static const unsigned MINMOD_LIMITER = 1;
  static const unsigned SUPERBEE_LIMITER = 2;
  static const unsigned ZERO_LIMITER = 3;
  static const unsigned EPSILON_LIMITER = 4;

// indices into primitive variables array
  static const unsigned RHO = 0;
  static const unsigned UX = 1;
  static const unsigned PR = 2;

// indicies into conserved variables array
  static const unsigned MX = 1;
  static const unsigned ER = 2;

/**
 * Minmod function for two parameters.
 */
  double minmod(double a1, double a2)
  {
    if (a1>0 && a2>0)
      return std::min(a1, a2);
    if (a1<0 && a2<0)
      return std::max(a1, a2);
    return 0.0;
  }

/**
 * Minmod function for two parameters.
 */
  double minmod(double a1, double a2, double a3)
  {
    if (a1>0 && a2>0 && a3>0)
      return Lucee::min3(a1, a2, a3);
    if (a1<0 && a2<0 && a3<0)
      return Lucee::max3(a1, a2, a3);
    return 0.0;
  }

/** Class id: this is used by registration system */
  const char *MusclHancock1DUpdater::id = "MusclHancock1D";

  MusclHancock1DUpdater::MusclHancock1DUpdater()
    : Lucee::UpdaterIfc() 
  {
  }

  void
  MusclHancock1DUpdater::declareTypes() 
  {
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }

  void
  MusclHancock1DUpdater::readInput(Lucee::LuaTable& tbl) 
  {
// call base class method
    UpdaterIfc::readInput(tbl);
// gas adiabatic constant
    gas_gamma = tbl.getNumber("gas_gamma");
// cfl number
    cfl = tbl.getNumber("cfl");

    limiter = AVERAGE_LIMITER;
    if (tbl.hasString("limiter"))
    {
      std::string lim = tbl.getString("limiter");
      if (lim == "average")
        limiter = AVERAGE_LIMITER;
      else if (lim == "minmod")
        limiter = MINMOD_LIMITER;
      else if (lim == "superbee")
        limiter = SUPERBEE_LIMITER;
      else if (lim == "zero")
        limiter = ZERO_LIMITER;
      else if (lim == "epsilon")
        limiter = EPSILON_LIMITER;
      else
      {
        Lucee::Except lce("MusclHancock1DUpdater::readInput: Do not recognize limiter type '");
        lce << lim << "'" << std::endl;
        throw lce;
      }
    }
  }

  void
  MusclHancock1DUpdater::initialize()
  {
// call base class initialization
    UpdaterIfc::initialize();

// get hold of grid
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();
// compute factor for use in EPSILON limiter
    double dx = grid.getDx(0);
    epsFac = dx*dx*dx;

// local region to index
    Lucee::Region<1, int> localRgn = grid.getLocalBox();

// ghost cells along slice
    int lg[1], ug[1];
    lg[0] = 2; ug[0] = 2;
// create region for allocating fields
    int lower[1], upper[1];
    lower[0] = localRgn.getLower(0); 
    upper[0] = localRgn.getUpper(0);
    Lucee::Region<1, int> slice(lower, upper);
// allocate memory (number of equations is 3)
    slopes = Lucee::Field<1, double>(slice, 3, lg, ug); // slopes
    predict = Lucee::Field<1, double>(slice, 3, lg, ug); // predicted solution
    prim = Lucee::Field<1, double>(slice, 3, lg, ug); // primitive variables
  }

  Lucee::UpdaterStatus
  MusclHancock1DUpdater::update(double t) 
  {
// get hold of grid
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();
// get input/output arrays
    const Lucee::Field<1, double>& q = this->getInp<Lucee::Field<1, double> >(0);
    Lucee::Field<1, double>& qNew = this->getOut<Lucee::Field<1, double> >(0);

// time-step
    double dt = t-this->getCurrTime();
    double dtdx = dt/grid.getDx(0);
// local region to index
    Lucee::Region<1, int> localRgn = grid.getLocalBox();

// pointers to data
    Lucee::ConstFieldPtr<double> pPtr = prim.createConstPtr();
    Lucee::ConstFieldPtr<double> plPtr = prim.createConstPtr();
    Lucee::ConstFieldPtr<double> prPtr = prim.createConstPtr();
    Lucee::FieldPtr<double> slpPtr = slopes.createPtr();

// indices of region to update
    int lower = localRgn.getLower(0);
    int upper = localRgn.getUpper(0);

    std::vector<double> ldiff(3), rdiff(3); // for jump in averages
    std::vector<double> ldelta(3), rdelta(3); // for projected jumps
    std::vector<double> projSlopes(3), testRecon(3);

// compute primitive variables from conserved variables
    calcPrimVars(q, prim);

    double emax = 0.0;
// compute maximum eigenvalue to check CFL condition
    for (int i=lower; i<upper; ++i)
    {
// attach pointer to primitive variables
      prim.setPtr(pPtr, i);
// compute maximum eigenvalue
      emax = std::max(emax, std::fabs(pPtr[UX]) + std::sqrt(gas_gamma*pPtr[PR]/pPtr[RHO]));
    }
// check if CFL condition is satisfied
    double cfla = dtdx*emax;
    if (cfla > 1.01*cfl) // the 1.01 fudge is required to avoid grinding away with same dt
      return Lucee::UpdaterStatus(false, dt*cfl/cfla);

// compute slopes for linear reconstruction in each cell (we need to
// compute slopes in one extra cell on either side so they can be used
// to update the first and last cell in the domain)
    for (int i=lower-1; i<upper+1; ++i)
    {
// attach pointers to left, current and right cells
      prim.setPtr(plPtr, i-1);
      prim.setPtr(pPtr, i);
      prim.setPtr(prPtr, i+1);

// compute jump in averages
      for (unsigned k=0; k<3; ++k)
      {
        ldiff[k] = pPtr[k]-plPtr[k];
        rdiff[k] = prPtr[k]-pPtr[k];
      }
// project these jumps onto left eigenvectors
      projectOnLeftEigenvectors(&pPtr[0], &ldiff[0], &ldelta[0]);
      projectOnLeftEigenvectors(&pPtr[0], &rdiff[0], &rdelta[0]);

// calculate slopes (these do not have the 1/dx term as it is taken
// into account in predictor step)
      for (unsigned k=0; k<3; ++k)
        projSlopes[k] = limaverage(ldelta[k], rdelta[k]);
// reconstruct slopes
      slopes.setPtr(slpPtr, i);
      reconWihRightEigenvectors(&pPtr[0], &projSlopes[0], &slpPtr[0]);
    }

// attach pointers for use in predictor step
    Lucee::ConstFieldPtr<double> cslpPtr = slopes.createConstPtr();
    Lucee::FieldPtr<double> prdPtr = predict.createPtr();

    double dtdx2 = 0.5*dtdx;
// compute predicted values of solution
    for (int i=lower-1; i<upper+1; ++i)
    {
      prim.setPtr(pPtr, i);
      slopes.setPtr(cslpPtr, i);
      predict.setPtr(prdPtr, i);

// compute predicted primitive variables at dt/2
      prdPtr[RHO] = pPtr[RHO] - dtdx2*(pPtr[RHO]*cslpPtr[UX] + pPtr[UX]*cslpPtr[RHO]); // density
      prdPtr[UX] = pPtr[UX] - dtdx2*(pPtr[UX]*cslpPtr[UX] + 1/pPtr[RHO]*cslpPtr[PR]); // velocity
      prdPtr[PR] = pPtr[PR] - dtdx2*(pPtr[UX]*cslpPtr[PR] + gas_gamma*pPtr[PR]*cslpPtr[UX]); // pressure
    }

// attach pointers for use in corrector step
    Lucee::ConstFieldPtr<double> cslplPtr = slopes.createConstPtr();
    Lucee::ConstFieldPtr<double> cslprPtr = slopes.createConstPtr();
    Lucee::FieldPtr<double> qNewlPtr = qNew.createPtr();
    Lucee::FieldPtr<double> qNewrPtr = qNew.createPtr();

    std::vector<double> ledge(3), redge(3); // primitive variables at left/right of edge
    std::vector<double> numFlux(3); // numerical flux

// qnew <- qold
    qNew.copy(q);
// compute corrected solution (this loop is over edges, hence we need
// one extra upper index)
    for (int i=lower; i<upper+1; ++i)
    {
      predict.setPtr(plPtr, i-1); // cell left of edge
      predict.setPtr(prPtr, i); // cell right of edge

      slopes.setPtr(cslplPtr, i-1); // cell left of edge
      slopes.setPtr(cslprPtr, i); // cell right of edge

// compute predicted solution at dt/2 at left/right of edge
      for (unsigned k=0; k<3; ++k)
      {
        ledge[k] = plPtr[k] + 0.5*cslplPtr[k]; // prim vars on left of edge
        redge[k] = prPtr[k] - 0.5*cslprPtr[k]; // prim vars on right of edge
      }

// evaluate numerical flux function using left/right states
      calcNumericalFlux(ledge, redge, numFlux);

      qNew.setPtr(qNewlPtr, i-1); // cell left of edge
      qNew.setPtr(qNewrPtr, i); // cell right of edge
// update solution in each cell connected to edge
      for (unsigned k=0; k<3; ++k)
      {
        qNewlPtr[k] += -dtdx*numFlux[k];
        qNewrPtr[k] += dtdx*numFlux[k];
      }
    }

    return Lucee::UpdaterStatus(true, dt*cfl/cfla);
  }

  double
  MusclHancock1DUpdater::limaverage(double a, double b)
  {
    double av = 0;
    switch (limiter)
    {
      case AVERAGE_LIMITER:
        av = 0.5*(a+b);
        break;
        
      case MINMOD_LIMITER:
        if (a*b > 0)
          av = minmod(0.5*(a+b), 2*a, 2*b);
        else
          av = 0.0;
        break;
          
      case SUPERBEE_LIMITER:
        av = 0.0;
        break;

      case ZERO_LIMITER:
        av = 0.0;
        break;

        case EPSILON_LIMITER:
        {
          double a2 = a*a, b2 = b*b, e2 = epsFac;
          av = ((b2+e2)*a + (a2+e2)*b)/(a2+b2+2*e2);
        }
        break;
          
        default:
// this can not happen
          break;
    };
    return av;
  }

  void
  MusclHancock1DUpdater::calcPrimVars(const Lucee::Field<1, double>& cv, Lucee::Field<1, double>& pv)
  {
    Lucee::ConstFieldPtr<double> cvPtr = cv.createConstPtr();
    Lucee::FieldPtr<double> pvPtr = pv.createPtr();

    for (int i=cv.getLowerExt(0); i<cv.getUpperExt(0); ++i)
    {
      cv.setPtr(cvPtr, i);
      pv.setPtr(pvPtr, i);

      pvPtr[RHO] = cvPtr[RHO]; // density
      pvPtr[UX] = cvPtr[MX]/cvPtr[RHO]; // velocity
      pvPtr[PR] = (cvPtr[ER] - 0.5*pvPtr[RHO]*pvPtr[UX]*pvPtr[UX])*(gas_gamma-1); // pressure
    }
  }

  void
  MusclHancock1DUpdater::calcNumericalFlux(const std::vector<double>& pvl, const std::vector<double> &pvr,
    std::vector<double>& nf)
  {
// compute sound speeds with left and right states
    double cdl = std::sqrt(gas_gamma*pvl[PR]/pvl[RHO]);
    double cdr = std::sqrt(gas_gamma*pvr[PR]/pvr[RHO]);
// compute speed needed to evaluate Lax fluxes
    double c = std::max(std::abs(pvl[UX])+cdl, std::abs(pvr[UX])+cdr);

    std::vector<double> fl(3), fr(3), cvl(3), cvr(3);
// calculate fluxes and conserved variables on left/right of edge
    calcFlux(pvl, fl); // left flux
    calcFlux(pvr, fr); // right flux
// (THIS IS NOT NEEDED AS WE KNOW CONSERVED VARS. SO PASS IN CONSERVED VARS AS PARAMETERS TO FUNCTION)
    calcConsVars(pvl, cvl); // left conserved variables
    calcConsVars(pvr, cvr); // right conserved variables

// compute Lax fluxes
    for (unsigned i=0; i<3; ++i)
      nf[i] = 0.5*(fl[i] + fr[i]) - 0.5*c*(cvr[i]-cvl[i]);
  }

  void
  MusclHancock1DUpdater::calcFlux(const std::vector<double>& pv, std::vector<double>& flux)
  {
    flux[0] = pv[RHO]*pv[UX]; // density flux
    flux[1] = pv[RHO]*pv[UX]*pv[UX] + pv[PR]; // momentum flux
    double E = pv[PR]/(gas_gamma-1) + 0.5*pv[RHO]*pv[UX]*pv[UX];
    flux[2] = pv[UX]*(E+pv[PR]); // energy flux
  }

  void
  MusclHancock1DUpdater::calcConsVars(const std::vector<double>& pv, std::vector<double>& cv)
  {
    cv[0] = pv[RHO]; // density
    cv[1] = pv[RHO]*pv[UX]; // momentum
    cv[2] = pv[PR]/(gas_gamma-1) + 0.5*pv[RHO]*pv[UX]*pv[UX]; // energy
  }

  void
  MusclHancock1DUpdater::projectOnLeftEigenvectors(const double *pv, const double *vec, double *coeff)
  {
// compute some needed quantities first
    double b = gas_gamma-1;
    double c = std::sqrt(gas_gamma*pv[PR]/pv[RHO]); // sound speed
    double u = pv[UX];
    double q2 = u*u;
    double th = 0.5*q2;
    double h = c*c/(gas_gamma-1) + 0.5*q2; // specific enthalpy
    double b2c2 = b/(2*c*c);

    coeff[0] = (th+u*c/b)*vec[0] + (-u-c/b)*vec[1] + vec[2];
    coeff[1] = (2*h-2*q2)*vec[0] + 2*u*vec[1] - 2*vec[2];
    coeff[2] = (th-u*c/b)*vec[0] + (-u+c/b)*vec[1] + vec[2];

// scale all of these
    for (unsigned i=0; i<3; ++i)
      coeff[i] = b2c2*coeff[i];
  }

  void
  MusclHancock1DUpdater::reconWihRightEigenvectors(const double *pv, const double *coeff, double *vec)
  {
// compute some needed quantities first
    double b = gas_gamma-1;
    double c = std::sqrt(gas_gamma*pv[PR]/pv[RHO]); // sound speed
    double u = pv[UX];
    double q2 = u*u;
    double th = 0.5*q2;
    double h = c*c/(gas_gamma-1) + 0.5*q2; // specific enthalpy
    double b2c2 = b/(2*c*c);

    vec[0] = coeff[0] + coeff[1] + coeff[2];
    vec[1] = (u-c)*coeff[0] + u*coeff[1] + (u+c)*coeff[2];
    vec[2] = (h-u*c)*coeff[0] + 0.5*q2*coeff[1] + (h+u*c)*coeff[2];
  }
}
