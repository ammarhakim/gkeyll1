/**
 * @file	LcSOLMaxwellianParameterCalc.cpp
 *
 * @brief	Extremely basic updater that determines next parameters to use as inputs for SOLMaxwellianAtNodeCalc
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLMaxwellianParameterCalc.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  const char *SOLMaxwellianParameterCalc::id = "SOLMaxwellianParameterCalc";

  SOLMaxwellianParameterCalc::SOLMaxwellianParameterCalc()
  {
  }

  void
  SOLMaxwellianParameterCalc::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("SOLMaxwellianParameterCalc::readInput: Must specify element to use using 'basis3d'");
  }

  void
  SOLMaxwellianParameterCalc::initialize()
  {
    UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus
  SOLMaxwellianParameterCalc::update(double t)
  {
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    const Lucee::Field<3, double>& mom1dir3TargetIn = this->getInp<Lucee::Field<3, double> >(0);
    const Lucee::Field<3, double>& temperatureTargetIn = this->getInp<Lucee::Field<3, double> >(1);
    const Lucee::Field<3, double>& mom1dir3LastGuessIn = this->getInp<Lucee::Field<3, double> >(2);
    const Lucee::Field<3, double>& temperatureLastGuessIn = this->getInp<Lucee::Field<3, double> >(3);
    const Lucee::Field<3, double>& mom1dir3NumericalIn = this->getInp<Lucee::Field<3, double> >(4);
    const Lucee::Field<3, double>& temperatureNumericalIn = this->getInp<Lucee::Field<3, double> >(5);
    // Output: Next guess inputs for SOLMaxwellianAtNodeCalc
    Lucee::Field<3, double>& mom1dir3Out = this->getOut<Lucee::Field<3, double> >(0);
    Lucee::Field<3, double>& temperatureOut = this->getOut<Lucee::Field<3, double> >(1);

    Lucee::ConstFieldPtr<double> mom1dir3TargetPtr = mom1dir3TargetIn.createConstPtr();
    Lucee::ConstFieldPtr<double> temperatureTargetPtr = temperatureTargetIn.createConstPtr();

    Lucee::ConstFieldPtr<double> mom1dir3LastGuessPtr = mom1dir3LastGuessIn.createConstPtr();
    Lucee::ConstFieldPtr<double> temperatureLastGuessPtr = temperatureLastGuessIn.createConstPtr();

    Lucee::ConstFieldPtr<double> mom1dir3NumericalPtr = mom1dir3NumericalIn.createConstPtr();
    Lucee::ConstFieldPtr<double> temperatureNumericalPtr = temperatureNumericalIn.createConstPtr();
    
    // Output to be used in next maxwellian calculation
    Lucee::FieldPtr<double> mom1dir3Ptr = mom1dir3Out.createPtr();
    Lucee::FieldPtr<double> temperaturePtr = temperatureOut.createPtr();

    unsigned nlocal3d = nodalBasis3d->getNumNodes();

    int idx[3];
    Lucee::Region<3, int> localRgn = grid.getLocalRegion();
    Lucee::RowMajorSequencer<3> seq(localRgn);

    const double tol = 1e-4;
    bool convergenceStatus = true;

    // Loop over local region
    while (seq.step())
    {
      seq.fillWithIndex(idx);

      // Set input fields to right location
      mom1dir3TargetIn.setPtr(mom1dir3TargetPtr, idx);
      temperatureTargetIn.setPtr(temperatureTargetPtr, idx);
      mom1dir3LastGuessIn.setPtr(mom1dir3LastGuessPtr, idx);
      temperatureLastGuessIn.setPtr(temperatureLastGuessPtr, idx);
      mom1dir3NumericalIn.setPtr(mom1dir3NumericalPtr, idx);
      temperatureNumericalIn.setPtr(temperatureNumericalPtr, idx);
      // Set output fields to right location
      mom1dir3Out.setPtr(mom1dir3Ptr, idx);
      temperatureOut.setPtr(temperaturePtr, idx);

      for (int i = 0; i < nlocal3d; i++)
      {
        mom1dir3Ptr[i] = (mom1dir3TargetPtr[i]/mom1dir3NumericalPtr[i])*mom1dir3LastGuessPtr[i];
        temperaturePtr[i] = (temperatureTargetPtr[i]/temperatureNumericalPtr[i])*temperatureLastGuessPtr[i];

        // Check convergence
        if (std::fabs(mom1dir3NumericalPtr[i] - mom1dir3TargetPtr[i]) > mom1dir3TargetPtr[i]*tol)
          convergenceStatus = false;
        
        if (std::fabs(temperatureNumericalPtr[i] - temperatureTargetPtr[i]) > temperatureTargetPtr[i]*tol)
          convergenceStatus = false;
      }
    }

    return Lucee::UpdaterStatus(convergenceStatus, 0.0);
  }

  void
  SOLMaxwellianParameterCalc::declareTypes()
  {
    // desired <v> = u*n
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // desired temperature (in joules)
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // last guess <v> = u*n
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // last guess temperature (in joules)
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // numerical <v> = u*n
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // numerical temperature (in joules)
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // next guess <v> = u*n
    this->appendOutVarType(typeid(Lucee::Field<5, double>));
    // next guess temperature (in joules)
    this->appendOutVarType(typeid(Lucee::Field<5, double>));
  }
}
