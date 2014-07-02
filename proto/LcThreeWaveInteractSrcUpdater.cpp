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
// fetch fields: in order a, b, f
    Lucee::Field<1, double>& a = this->getOut<Lucee::Field<1, double> >(0);
    Lucee::Field<1, double>& b = this->getOut<Lucee::Field<1, double> >(1);
    Lucee::Field<1, double>& f = this->getOut<Lucee::Field<1, double> >(2);

    Lucee::FieldPtr<double> aPtr = a.createPtr();
    Lucee::FieldPtr<double> bPtr = b.createPtr();
    Lucee::FieldPtr<double> fPtr = f.createPtr();

    unsigned maxIter = 0; // for diagnostics

    int idx[1];
// local region to index
    Lucee::Region<1, int> localRgn = grid.getLocalRegion();
    Lucee::RowMajorSequencer<1> seq(localRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      int ix = idx[0];
      a.setPtr(aPtr, ix);
      b.setPtr(bPtr, ix);
      f.setPtr(fPtr, ix);

      std::complex<double> ac(aPtr[0], aPtr[1]);
      std::complex<double> bc(bPtr[0], bPtr[1]);
      std::complex<double> fc(fPtr[0], fPtr[1]);

      std::complex<double> acp, bcp, fcp; // updated values

      unsigned nIter = 0;
// NOTE: This loop solves the ODE implicitly using an iterative
// scheme. A backward Euler scheme is used. This may diffusive.
// Hence, it may be better to use the mid-point rule or even solve the
// system using Boost ODE solvers. Will change if needed.
      bool done = false;
      while (!done)
      {
        nIter++;
// first compute predicted value of 'f'
        fcp = fc - dt*ac*std::conj(bc);
// update 'a' and 'b' using an implicit update (one can invert system
// analytically) using predicted value of f
        double fc1 = 1/(1+dt*dt*std::norm(fcp));
        acp = (ac-dt*bc*fcp)*fc1;
        bcp = (bc+dt*ac*std::conj(fcp))*fc1;
// copy back before loop
        ac = acp; bc = bcp; fc = fcp;

// check for convergence
        std::complex<double> fcpp = fc - dt*ac*std::conj(bc);
        if (epsCmp(fcp.real(), fcpp.real(), relTol))
          done = true;
      }
      maxIter = nIter > maxIter ? nIter : maxIter;

// copy updated values into appropriate fields
      aPtr[0] = acp.real(); aPtr[1] = acp.imag();
      bPtr[0] = bcp.real(); bPtr[1] = bcp.imag();
      fPtr[0] = fcp.real(); fPtr[1] = fcp.imag();
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
}

