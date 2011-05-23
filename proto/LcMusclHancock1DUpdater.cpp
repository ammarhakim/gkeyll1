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

namespace Lucee
{
  static const unsigned AVERAGE_LIMITER = 0;
  static const unsigned MINMOD_LIMITER = 1;
  static const unsigned SUPERBEE_LIMITER = 2;
  static const unsigned ZERO_LIMITER = 3;

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
    slopes = Lucee::Field<1, double>(slice, 3, lg, ug);
    predict = Lucee::Field<1, double>(slice, 3, lg, ug);

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

// pointers to data
    Lucee::ConstFieldPtr<double> qPtr = q.createConstPtr();
    Lucee::ConstFieldPtr<double> qlPtr = q.createConstPtr();
    Lucee::ConstFieldPtr<double> qrPtr = q.createConstPtr();
    Lucee::FieldPtr<double> slpPtr = slopes.createPtr();

// indices of region to update
    int lower = localRgn.getLower(0);
    int upper = localRgn.getUpper(0);

// arrays to store primitive variables
    std::vector<double> pvl(3), pv(3), pvr(3);

// compute slopes for linear reconstruction in each cell
    for (unsigned i=lower; i<upper; ++i)
    {
// attach pointers to left, current and right cells (could be done
// more efficiently as the primitive variables in this implementation
// are being computed too many times)
      q.setPtr(qlPtr, i-1);
      q.setPtr(qPtr, i);
      q.setPtr(qrPtr, i+1);

// compute primitive variables in each cell
      calcPrimVars(&qlPtr[0], &pvl[0]);
      calcPrimVars(&qPtr[0], &pv[0]);
      calcPrimVars(&qrPtr[0], &pvr[0]);

// calculate slope (these do not have the 1/dx term as it is taken
// into account in predictor step)
      slopes.setPtr(slpPtr, i);
      for (unsigned i=0; i<3; ++i)
        slpPtr[i] = limaverage(pv[i]-pvl[i], pvr[i]-pv[i]);
    }

// attach pointers for use in predictor stop
    Lucee::ConstFieldPtr<double> cslpPtr = slopes.createConstPtr();
    Lucee::FieldPtr<double> prdPtr = predict.createPtr();

    double dtdx = dt/grid.getDx(0);
    double dtdx2 = 0.5*dtdx;
// compute predicted values of solution
    for (unsigned i=lower; i<upper; ++i)
    {
      q.setPtr(qPtr, i);
      slopes.setPtr(cslpPtr, i);
      predict.setPtr(prdPtr, i);

// compute primitive variables first
      calcPrimVars(&qPtr[0], &pv[0]);

// compute predicted primitive variables at dt/2
      prdPtr[0] = pv[0] - dtdx2*(pv[0]*cslpPtr[1] + pv[1]*cslpPtr[0]); // density
      prdPtr[1] = pv[1] - dtdx2*(pv[1]*cslpPtr[1] + 1/pv[0]*cslpPtr[2]); // velocity
      prdPtr[2] = pv[2] - dtdx2*(pv[1]*cslpPtr[2] + gas_gamma*pv[2]*cslpPtr[1]); // pressure
    }

// attach pointers for use in corrector step

// compute corrected solution (this loop is over edges, hence we need
// one extra upper index)
    for (unsigned i=lower; i<upper+1; ++i)
    {
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
    switch (limiter)
    {
        case AVERAGE_LIMITER:
          return 0.5*(a+b);
          break;

        case MINMOD_LIMITER:
          break;
          
        case SUPERBEE_LIMITER:
          break;

        case ZERO_LIMITER:
          return 0;
          break;

        default:
// this can not happen
          break;
    };
    return 0;
  }

  void
  MusclHancock1DUpdater::calcPrimVars(const double *cv, double *pv)
  {
    pv[0] = cv[0]; // density
    pv[1] = cv[1]/cv[0]; // velocity
    pv[2] = (cv[2] - 0.5*cv[1]*cv[1]/cv[0])*(gas_gamma-1);
  }
}
