/**
 * @file	LcThreeWaveInteractModSrcUpdater.cpp
 *
 * @brief	Three wave interaction source updater.
 */

// lucee includes
#include <LcThreeWaveInteractModSrcUpdater.h>
#include <LcStructGridField.h>

// eigen inlcudes
#include <Eigen/Eigen>

// boost includes
#include <boost/numeric/odeint.hpp>
//using namespace boost::numeric::odeint;
// std includes
#include <complex>

namespace Lucee
{
// set ids for module system
  const char *ThreeWaveInteractModSrcUpdater::id = "ThreeWaveInteractModSrc";

  bool
  ThreeWaveInteractModSrcUpdater::epsCmp(double a, double b) const
  {
// The 1e3 seems arbitrary, but prevents a divide by zero
    if (fabs(a) <= 1e3*std::numeric_limits<double>::epsilon());
      return fabs(a-b) <= relTol;
    return fabs(1-b/a) <= relTol;
  }

  void
  ThreeWaveInteractModSrcUpdater::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    std::vector<double> clst = tbl.getNumVec("constants");
    if (clst.size() != 3)
    {
      Lucee::Except lce(
        "ThreeWaveInteractModSrcUpdater::readInput: 'constants' should have exactly 3 entries. ");
      lce << clst.size() << " provided instead";
      throw lce;
    }

// for now assume constants are real (should be easy to fix if needed)
    c[0] = std::complex<double>(clst[0], 0.0);
    c[1] = std::complex<double>(clst[1], 0.0);
    c[2] = std::complex<double>(clst[2], 0.0);

    intType = TWI_IMPLICIT; // default integrator
// get integrator type (if specified)
    if (tbl.getString("integrator") == "rk4")
      intType = TWI_RK4;
    else if (tbl.getString("integrator") == "implicit")
      intType = TWI_IMPLICIT;

    relTol = 1.0e-4;
    if (intType == TWI_IMPLICIT)
      relTol = tbl.getNumber("relativeTol");
  }

  void
  ThreeWaveInteractModSrcUpdater::initialize()
  {
    UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus
  ThreeWaveInteractModSrcUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();

    double dt = t-this->getCurrTime();
// fetch fields: in order E1, E2, E3, W
    Lucee::Field<1, double>& e1 = this->getOut<Lucee::Field<1, double> >(0);
    Lucee::Field<1, double>& e2 = this->getOut<Lucee::Field<1, double> >(1);
    Lucee::Field<1, double>& e3 = this->getOut<Lucee::Field<1, double> >(2);
    Lucee::Field<1, double>& w = this->getOut<Lucee::Field<1, double> >(3);

    Lucee::FieldPtr<double> e1Ptr = e1.createPtr();
    Lucee::FieldPtr<double> e2Ptr = e2.createPtr();
    Lucee::FieldPtr<double> e3Ptr = e3.createPtr();
    Lucee::FieldPtr<double> wPtr = w.createPtr();

    std::vector<std::complex<double> > inpE(4);
    std::vector<std::complex<double> > outE(4);

    unsigned maxIter = 0; // for diagnostics

    int idx[1];
// local region to index
    Lucee::Region<1, int> localRgn = grid.getLocalRegion();
    Lucee::RowMajorSequencer<1> seq(localRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      int ix = idx[0];
      e1.setPtr(e1Ptr, ix);
      e2.setPtr(e2Ptr, ix);
      e3.setPtr(e3Ptr, ix);
      w.setPtr(wPtr, ix);

      inpE[0] = std::complex<double> (e1Ptr[0], e1Ptr[1]);
      inpE[1] = std::complex<double> (e2Ptr[0], e2Ptr[1]);
      inpE[2] = std::complex<double> (e3Ptr[0], e3Ptr[1]);
      inpE[3] = std::complex<double> (wPtr[0], wPtr[1]);

      unsigned nIter = 1;
      if (intType == TWI_IMPLICIT)
        nIter = stepImplicit(dt, inpE, outE);
      else if (intType == TWI_RK4)
        nIter = stepRK4(dt, inpE, outE);
      else
      { /** Can't happen */ }

// store maximum number of iterations as diagnostic
      maxIter = nIter > maxIter ? nIter : maxIter;

// copy updated values into appropriate fields
      e1Ptr[0] = outE[0].real(); e1Ptr[1] = outE[0].imag();
      e2Ptr[0] = outE[1].real(); e2Ptr[1] = outE[1].imag();
      e3Ptr[0] = outE[2].real(); e3Ptr[1] = outE[2].imag();
      wPtr[0]  = outE[3].real(); wPtr[1]  = outE[3].imag();
    }

    // if (maxIter>1)
    //   std::cout << "Max number of iterations " << maxIter << std::endl;
    
    return Lucee::UpdaterStatus();
  }
  

  void
  ThreeWaveInteractModSrcUpdater::declareTypes()
  {
// outputs a, b, f which are assumed to be complex fields
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }

  unsigned
  ThreeWaveInteractModSrcUpdater::stepImplicit(double dt, const std::vector<std::complex<double> >& inp,
    std::vector<std::complex<double> >& out)
  {
    throw Lucee::Except("Lucee::ThreeWaveInteractModSrcUpdater: Implicit update not implemented!");
    return 0;
  }

  unsigned
  ThreeWaveInteractModSrcUpdater::stepRK4(double dt, const std::vector<std::complex<double> >& inp,
    std::vector<std::complex<double> >& out)
  {
    boost::numeric::odeint::runge_kutta4<twstate> stepper;
    stepper.do_step( *this, inp, 0, out, dt);
    return 1;
  }

  void
  ThreeWaveInteractModSrcUpdater::operator() (const twstate &x, twstate &dxdt, const double t)
  {
    dxdt[0] = c[0]*std::conj(x[1])*std::conj(x[2]);
    dxdt[1] = c[1]*std::conj(x[0])*std::conj(x[2]);
    dxdt[3] = x[3];
    dxdt[4] = c[2]*std::conj(x[0])*std::conj(x[1]);
  }
}

