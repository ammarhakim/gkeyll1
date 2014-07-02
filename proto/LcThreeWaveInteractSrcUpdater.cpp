/**
 * @file	LcThreeWaveInteractSrcUpdater.cpp
 *
 * @brief	Implicit updater for 5-moment source terms
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

// set ids for module system
  const char *ThreeWaveInteractSrcUpdater::id = "ThreeWaveInteractSrc";

  void
  ThreeWaveInteractSrcUpdater::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);
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

// copy updated values into appropriate fields
      aPtr[0] = acp.real(); aPtr[1] = acp.imag();
      bPtr[0] = bcp.real(); bPtr[1] = bcp.imag();
      fPtr[0] = fcp.real(); fPtr[1] = fcp.imag();
    }
    
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

