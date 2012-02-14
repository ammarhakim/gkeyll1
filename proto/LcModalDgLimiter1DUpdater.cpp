/**
 * @file	LcModalDgLimiter1DUpdater.cpp
 *
 * @brief	Updater to apply limiters to DG solution.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAlignedRectCoordSys.h>
#include <LcMathLib.h>
#include <LcModalDgLimiter1DUpdater.h>
#include <LcStructuredGridBase.h>

namespace Lucee
{
/** Class id: this is used by registration system */
  const char *ModalDgLimiter1DUpdater::id = "ModalDgLimiter1D";

/**
 * Minmod function for three parameters.
 */
  static double minmod(double a1, double a2, double a3)
  {
    if (a1>0 && a2>0 && a3>0)
      return Lucee::min3(a1, a2, a3);
    if (a1<0 && a2<0 && a3<0)
      return Lucee::max3(a1, a2, a3);
    return 0.0;
  }

  double
  ModalDgLimiter1DUpdater::modifiedMinMod(double a, double b, double c, double dx) const
  {
    if (std::fabs(a) < Mfact*dx*dx)
      return a;
    else
      return minmod(a, b, c);
  }

  ModalDgLimiter1DUpdater::ModalDgLimiter1DUpdater()
    : UpdaterIfc(), numBasis(1)
  {
  }

  void
  ModalDgLimiter1DUpdater::readInput(Lucee::LuaTable& tbl)
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

// factor for slope correction
    Mfact = 0.0;
    if (tbl.hasNumber("M"))
      Mfact = tbl.getNumber("M");
  }

  void
  ModalDgLimiter1DUpdater::initialize()
  {
// call base class method
    Lucee::UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus
  ModalDgLimiter1DUpdater::update(double t)
  {
    if (numBasis < 2)
// no need to apply limiters to first-order scheme
      return Lucee::UpdaterStatus();

// get hold of grid
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();
// get output arrays
    Lucee::Field<1, double>& q = this->getOut<Lucee::Field<1, double> >(0);

// local region to index
    Lucee::Region<1, int> localRgn = grid.getLocalRegion();

// create coordinate system along X-axis
    Lucee::AlignedRectCoordSys coordSys(0);

    if ( (q.getNumComponents() != meqn*numBasis) && (q.getNumComponents() != meqn*numBasis) )
    {
      Lucee::Except lce(
        "ModalDg1DUpdater::update: Number of components in input/output fields should be ");
      lce << meqn*numBasis << ". Instead provided " << q.getNumComponents() << std::endl;
    }

// lower and upper bounds of slice.
    int sliceLower = localRgn.getLower(0);
    int sliceUpper = localRgn.getUpper(0);

    double dx, xc[3];

// space to store differences and projections
    std::vector<double> fDiff(meqn), bDiff(meqn);
    std::vector<double> linProj(meqn), fProj(meqn), bProj(meqn);
    std::vector<double> limSlopes(meqn);

    Lucee::FieldPtr<double> qPtr = q.createPtr();
    Lucee::FieldPtr<double> qPtrP = q.createPtr();
    Lucee::FieldPtr<double> qPtrM = q.createPtr();

    std::vector<bool> changed(meqn); // flags to see if linear term was modified

    for (unsigned i=sliceLower; i<sliceUpper; ++i)
    {
      q.setPtr(qPtr, i); // cell i
      q.setPtr(qPtrP, i+1); // cell i+1
      q.setPtr(qPtrM, i-1); // cell i-1

// get cell spacing
      grid.setIndex(i);
      dx = grid.getDx(0);

// compute forward and backward differences of averages
      for (unsigned k=0; k<meqn; ++k)
      {
        fDiff[k] = qPtrP[k] - qPtr[k];
        bDiff[k] = qPtr[k] - qPtrM[k];
      }

// compute projection of linear and forward/backward terms: these
// calls essentially say to compute flux Jacobian using cell averages
// and using that for the projection. Note that the linear
// coefficients are stored after cell averages and hence the q[0+meqn]
// in the first call

      equation->projectOnLeftEigenvectors(coordSys, qPtr, &qPtr[0+meqn], &linProj[0]);
      equation->projectOnLeftEigenvectors(coordSys, qPtr, &fDiff[0], &fProj[0]);
      equation->projectOnLeftEigenvectors(coordSys, qPtr, &bDiff[0], &bProj[0]);

// now compute limited slopes
      for (unsigned k=0; k<meqn; ++k)
      {
        limSlopes[k] = modifiedMinMod(linProj[k], fProj[k], bProj[k], dx);
        if (limSlopes[k] == linProj[k])
          changed[k] = false;
        else
          changed[k] = true;
      }

// reconstruct slopes
      equation->reconWithRightEigenvectors(coordSys, qPtr, &limSlopes[0], &qPtr[0+meqn]);

// check if slope was modified and if it was, zap high-order
// coefficients
      if (numBasis > 2)
      {
        for (unsigned k=0; k<meqn; ++k)
        {
          if (changed[k])
            for (unsigned c=2; c<numBasis; ++c)
              qPtr[meqn*c+k] = 0.0;
        }
      }
    }

    return Lucee::UpdaterStatus();
  }

  void
  ModalDgLimiter1DUpdater::declareTypes()
  {
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }
}
