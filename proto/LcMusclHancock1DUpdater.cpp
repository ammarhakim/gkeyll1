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
#include <LcMusclHancock1DUpdater.h>
#include <LcStructuredGridBase.h>

// std includes
#include <cmath>

namespace Lucee
{
  static const unsigned AVERAGE_LIMITER = 0;
  static const unsigned MINMOD_LIMITER = 1;
  static const unsigned SUPERBEE_LIMITER = 2;
  static const unsigned ZERO_LIMITER = 3;

  static const unsigned RHO = 0;
  static const unsigned UX = 1;
  static const unsigned PR = 2;
  

/** Class id: this is used by registration system */
  const char *MusclHancock1DUpdater::id = "MusclHancock1D";

  MusclHancock1DUpdater::MusclHancock1DUpdater()
    : Lucee::UpdaterIfc() 
  {
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
      else if (lim == "min-mod")
        limiter = MINMOD_LIMITER;
      else if (lim == "superbee")
        limiter = SUPERBEE_LIMITER;
      else if (lim == "zero")
        limiter = ZERO_LIMITER;
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
// local region to index
    Lucee::Region<1, int> localRgn = grid.getLocalBox();

// compute primitive variables from conserved variables
    calcPrimVars(q, prim);

// pointers to data
    Lucee::ConstFieldPtr<double> pPtr = prim.createConstPtr();
    Lucee::ConstFieldPtr<double> plPtr = prim.createConstPtr();
    Lucee::ConstFieldPtr<double> prPtr = prim.createConstPtr();
    Lucee::FieldPtr<double> slpPtr = slopes.createPtr();

// indices of region to update
    int lower = localRgn.getLower(0);
    int upper = localRgn.getUpper(0);

// compute slopes for linear reconstruction in each cell
    for (int i=lower; i<upper; ++i)
    {
// attach pointers to left, current and right cells
      prim.setPtr(plPtr, i-1);
      prim.setPtr(pPtr, i);
      prim.setPtr(prPtr, i+1);

// calculate slope (these do not have the 1/dx term as it is taken
// into account in predictor step)
      slopes.setPtr(slpPtr, i);
      for (unsigned k=0; k<3; ++k)
        slpPtr[k] = limaverage(pPtr[k]-plPtr[k], prPtr[k]-pPtr[k]);
    }

// attach pointers for use in predictor step
    Lucee::ConstFieldPtr<double> cslpPtr = slopes.createConstPtr();
    Lucee::FieldPtr<double> prdPtr = predict.createPtr();

    double dtdx = dt/grid.getDx(0);
    double dtdx2 = 0.5*dtdx;
// compute predicted values of solution
    for (int i=lower; i<upper; ++i)
    {
      prim.setPtr(pPtr, i);
      slopes.setPtr(cslpPtr, i);
      predict.setPtr(prdPtr, i);

// compute predicted primitive variables at dt/2
      prdPtr[0] = pPtr[0] - dtdx2*(pPtr[0]*cslpPtr[1] + pPtr[1]*cslpPtr[0]); // density
      prdPtr[1] = pPtr[1] - dtdx2*(pPtr[1]*cslpPtr[1] + 1/pPtr[0]*cslpPtr[2]); // velocity
      prdPtr[2] = pPtr[2] - dtdx2*(pPtr[1]*cslpPtr[2] + gas_gamma*pPtr[2]*cslpPtr[1]); // pressure
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
      prim.setPtr(plPtr, i-1); // cell left of edge
      prim.setPtr(prPtr, i); // cell right of edge

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

// update solution in each cell connected to edge
      qNew.setPtr(qNewlPtr, i-1); // cell left of edge
      qNew.setPtr(qNewrPtr, i); // cell right of edge

      for (unsigned k=0; k<3; ++k)
      {
        qNewlPtr[k] += -dtdx*numFlux[k];
        qNewrPtr[k] += dtdx*numFlux[k];
      }
    }

    return Lucee::UpdaterStatus();
  }

  void
  MusclHancock1DUpdater::declareTypes() 
  {
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
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
          av = 0.0;
          break;
          
      case SUPERBEE_LIMITER:
          av = 0.0;
          break;

      case ZERO_LIMITER:
          av = 0.0;
          break;
          
        default:
// this can not happen
          break;
    };
    return av;
  }

  void
  MusclHancock1DUpdater::calcPrimVars(const Lucee::Field<1, double>& cv, Lucee::Field<1, double> &pv)
  {
    Lucee::ConstFieldPtr<double> cvPtr = cv.createConstPtr();
    Lucee::FieldPtr<double> pvPtr = pv.createPtr();

    for (int i=cv.getLowerExt(0); i<cv.getUpperExt(0); ++i)
    {
      cv.setPtr(cvPtr, i);
      pv.setPtr(pvPtr, i);

      pvPtr[0] = cvPtr[0]; // density
      pvPtr[1] = cvPtr[1]/cvPtr[0]; // velocity
      pvPtr[2] = (cvPtr[2] - 0.5*cvPtr[1]*cvPtr[1]/cvPtr[0])*(gas_gamma-1); // pressure
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
    calcConsVars(pvl, cvl); // left conserved variables
    calcConsVars(pvr, cvr); // right conserved variables

// compute Lax fluxes
    for (unsigned i=0; i<3; ++i)
      nf[i] = 0.5*(fl[i] + fr[i]) + 0.5*c*(cvr[i]-cvl[i]);
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
}
