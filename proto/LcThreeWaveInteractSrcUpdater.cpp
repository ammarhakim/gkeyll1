/**
 * @file	LcThreeWaveInteractSrcUpdater.cpp
 *
 * @brief	Three wave interaction source updater.
 */

// lucee includes
#include <LcThreeWaveInteractSrcUpdater.h>
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
  static 
  bool epsCmp(double a, double b, double epsDiff)
  {
// The 1e3 seems arbitrary, but prevents a divide by zero
    if (fabs(a) <= 1e3*std::numeric_limits<double>::epsilon());
      return fabs(a-b) <= epsDiff;
    return fabs(1-b/a) <= epsDiff;
  }

// set ids for module system
  const char *ThreeWaveInteractSrcUpdater::id = "ThreeWaveInteractSrc";

  void
  ThreeWaveInteractSrcUpdater::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    std::vector<double> clst = tbl.getNumVec("constants");
    if (clst.size() != 3)
    {
      Lucee::Except lce(
        "ThreeWaveInteractSrcUpdater::readInput: 'constants' should have exactly 3 entries. ");
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
  ThreeWaveInteractSrcUpdater::initialize()
  {
    UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus
  ThreeWaveInteractSrcUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();

    double dt = t-this->getCurrTime();
// fetch fields: in order E1, E2, E3
    Lucee::Field<1, double>& e1 = this->getOut<Lucee::Field<1, double> >(0);
    Lucee::Field<1, double>& e2 = this->getOut<Lucee::Field<1, double> >(1);
    Lucee::Field<1, double>& e3 = this->getOut<Lucee::Field<1, double> >(2);

    Lucee::FieldPtr<double> e1Ptr = e1.createPtr();
    Lucee::FieldPtr<double> e2Ptr = e2.createPtr();
    Lucee::FieldPtr<double> e3Ptr = e3.createPtr();

    std::vector<std::complex<double> > inpE(3);
    std::vector<std::complex<double> > outE(3);

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

      inpE[0] = std::complex<double> (e1Ptr[0], e1Ptr[1]);
      inpE[1] = std::complex<double> (e2Ptr[0], e2Ptr[1]);
      inpE[2] = std::complex<double> (e3Ptr[0], e3Ptr[1]);

      if (intType == TWI_IMPLICIT)
        stepImplicit(dt, inpE, outE);
      else if (intType == TWI_RK4)
        stepRK4(dt, inpE, outE);
      else
      { /** Can't happen */ }

// copy updated values into appropriate fields
      e1Ptr[0] = outE[0].real(); e1Ptr[1] = outE[0].imag();
      e2Ptr[0] = outE[1].real(); e2Ptr[1] = outE[1].imag();
      e3Ptr[0] = outE[2].real(); e3Ptr[1] = outE[2].imag();
    }

    if (maxIter>1)
      std::cout << "Max number of iterations " << maxIter << std::endl;
    
    return Lucee::UpdaterStatus();
  }
  

  void
  ThreeWaveInteractSrcUpdater::declareTypes()
  {
// outputs a, b, f which are assumed to be complex fields
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }

  unsigned
  ThreeWaveInteractSrcUpdater::stepImplicit(double dt, const std::vector<std::complex<double> >& inp,
    std::vector<std::complex<double> >& out)
  {
    unsigned nIter = 0;
// // NOTE: This loop solves the ODE implicitly using an iterative
// // scheme. A backward Euler scheme is used. This may diffusive.
// // Hence, it may be better to use the mid-point rule or even solve the
// // system using Boost ODE solvers. Will change if needed.
//     bool done = false;
//       while (!done)
//       {
//         nIter++;
// // first compute predicted value of 'f'
//         fcp = fc - dt*ac*std::conj(bc);
// // update 'a' and 'b' using an implicit update (one can invert system
// // analytically) using predicted value of f
//         double fc1 = 1/(1+dt*dt*std::norm(fcp));
//         acp = (ac-dt*bc*fcp)*fc1;
//         bcp = (bc+dt*ac*std::conj(fcp))*fc1;
// // copy back before loop
//         ac = acp; bc = bcp; fc = fcp;

// // check for convergence
//         std::complex<double> fcpp = fc - dt*ac*std::conj(bc);
//         if (epsCmp(fcp.real(), fcpp.real(), relTol))
//           done = true;
//       }
    return nIter;
  }

  unsigned
  ThreeWaveInteractSrcUpdater::stepRK4(double dt, const std::vector<std::complex<double> >& inp,
    std::vector<std::complex<double> >& out)
  {
    throw Lucee::Except("ThreeWaveInteractSrcUpdater::stepRK4: Not implemented!");
    return 1;
  }
}

