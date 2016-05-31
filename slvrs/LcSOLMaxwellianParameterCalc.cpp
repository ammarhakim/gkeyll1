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

    const double tol = 1e-3;
    bool convergenceStatus = true;

    // Loop over local region
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      //std::cout << "idx[" << idx[0] << "," << idx[1] << "," << idx[2] << "]" << std::endl;

      // Set input fields to right location
      mom1dir3TargetIn.setPtr(mom1dir3TargetPtr, idx);
      mom1dir3LastGuessIn.setPtr(mom1dir3LastGuessPtr, idx);
      mom1dir3NumericalIn.setPtr(mom1dir3NumericalPtr, idx);
      
      temperatureTargetIn.setPtr(temperatureTargetPtr, idx);
      temperatureLastGuessIn.setPtr(temperatureLastGuessPtr, idx);
      temperatureNumericalIn.setPtr(temperatureNumericalPtr, idx);
      // Set output fields to right location
      mom1dir3Out.setPtr(mom1dir3Ptr, idx);
      temperatureOut.setPtr(temperaturePtr, idx);

      for (int i = 0; i < nlocal3d; i++)
      {

        //if (temperaturePtr[i] == 0.0)
        //  temperaturePtr[i] = 0.0
        //else
        temperaturePtr[i] = (temperatureTargetPtr[i]-temperatureNumericalPtr[i]) + temperatureLastGuessPtr[i];

        if (mom1dir3NumericalPtr[i] == 0.0)
          mom1dir3Ptr[i] = 0.0;
        else
          mom1dir3Ptr[i] = (mom1dir3TargetPtr[i]-mom1dir3NumericalPtr[i]) + mom1dir3LastGuessPtr[i];

        /*if (std::isnan(mom1dir3NumericalPtr[i]) || std::isinf(mom1dir3NumericalPtr[i]))
        {
        std::cout << "mom1dir3TargetPtr " << i << " = " << mom1dir3TargetPtr[i] << std::endl;
        std::cout << "mom1dir3NumericalPtr " << i << " = " << mom1dir3NumericalPtr[i] << std::endl;
        std::cout << "mom1dir3LastGuessPtr " << i << " = " << mom1dir3LastGuessPtr[i] << std::endl << std::endl;
        }

        if (std::isnan(temperatureNumericalPtr[i]) || std::isinf(temperatureNumericalPtr[i]))
        {
        std::cout << "temperatureTargetPtr " << i << " = " << temperatureTargetPtr[i] << std::endl;
        std::cout << "temperatureNumericalPtr " << i << " = " << temperatureNumericalPtr[i] << std::endl;
        std::cout << "temperatureLastGuessPtr " << i << " = " << temperatureLastGuessPtr[i] << std::endl << std::endl;
        }
        */
        // Check convergence
        if (std::fabs(mom1dir3NumericalPtr[i] - mom1dir3TargetPtr[i]) > std::fabs(mom1dir3TargetPtr[i]*tol) )
        {
          //if (mom1dir3NumericalPtr[i] == 0.0)
          //  mom1dir3Ptr[i] = 0.0;
          //else
          //  mom1dir3Ptr[i] = (mom1dir3TargetPtr[i]-mom1dir3NumericalPtr[i]) + mom1dir3LastGuessPtr[i];
            //mom1dir3Ptr[i] = (mom1dir3TargetPtr[i]/mom1dir3NumericalPtr[i])*mom1dir3LastGuessPtr[i];

          //if (std::fabs(mom1dir3NumericalPtr[i] - mom1dir3TargetPtr[i]) > std::fabs(tol*mom1dir3TargetPtr[i]) )
          /*{
            std::cout << "idx[" << idx[0] << "," << idx[1] << "," << idx[2] << "]" << std::endl;
            std::cout << "mom1dir3TargetPtr " << i << " = " << mom1dir3TargetPtr[i] << std::endl;
            std::cout << "mom1dir3NumericalPtr " << i << " = " << mom1dir3NumericalPtr[i] << std::endl;
            std::cout << "mom1dir3LastGuessPtr " << i << " = " << mom1dir3LastGuessPtr[i] << std::endl;
            std::cout << "mom1dirNextGuess " << i << " = " << mom1dir3Ptr[i] << std::endl;
            std::cout << "temperatureTargetPtr " << i << " = " << temperatureTargetPtr[i] << std::endl;
            std::cout << "temperatureNumericalPtr " << i << " = " << temperatureNumericalPtr[i] << std::endl;
            std::cout << "temperatureLastGuessPtr " << i << " = " << temperatureLastGuessPtr[i] << std::endl;
            std::cout << "temperatureNextGuess " << i << " = " << temperaturePtr[i] << std::endl;
            std::cout << "relErrorU = " << std::fabs( (mom1dir3NumericalPtr[i] - mom1dir3TargetPtr[i])/mom1dir3TargetPtr[i] ) << std::endl;
            std::cout << "relErrorT = " << std::fabs( (temperatureNumericalPtr[i] - temperatureTargetPtr[i])/temperatureTargetPtr[i] ) << std::endl << std::endl;
          }*/
          convergenceStatus = false;
        }
        //else 
        //{
          // Don't modify guess if we are already converged
        //  mom1dir3Ptr[i] = mom1dir3LastGuessPtr[i];
        //}
        
        if (std::fabs(temperatureNumericalPtr[i] - temperatureTargetPtr[i]) > std::fabs(temperatureTargetPtr[i]*tol) )
        {
          //temperaturePtr[i] = (temperatureTargetPtr[i]/temperatureNumericalPtr[i])*temperatureLastGuessPtr[i];
          //temperaturePtr[i] = (temperatureTargetPtr[i]-temperatureNumericalPtr[i]) + temperatureLastGuessPtr[i];
          //if (std::fabs(mom1dir3NumericalPtr[i] - mom1dir3TargetPtr[i]) > std::fabs(tol*mom1dir3TargetPtr[i]) )
          /*{
            std::cout << "idx[" << idx[0] << "," << idx[1] << "," << idx[2] << "]" << std::endl;
            std::cout << "mom1dir3TargetPtr " << i << " = " << mom1dir3TargetPtr[i] << std::endl;
            std::cout << "mom1dir3NumericalPtr " << i << " = " << mom1dir3NumericalPtr[i] << std::endl;
            std::cout << "mom1dir3LastGuessPtr " << i << " = " << mom1dir3LastGuessPtr[i] << std::endl;
            std::cout << "mom1dirNextGuess " << i << " = " << mom1dir3Ptr[i] << std::endl;
            std::cout << "temperatureTargetPtr " << i << " = " << temperatureTargetPtr[i] << std::endl;
            std::cout << "temperatureNumericalPtr " << i << " = " << temperatureNumericalPtr[i] << std::endl;
            std::cout << "temperatureLastGuessPtr " << i << " = " << temperatureLastGuessPtr[i] << std::endl;
            std::cout << "temperatureNextGuess " << i << " = " << temperaturePtr[i] << std::endl;
            std::cout << "relErrorU = " << std::fabs( (mom1dir3NumericalPtr[i] - mom1dir3TargetPtr[i])/mom1dir3TargetPtr[i] ) << std::endl;
            std::cout << "relErrorT = " << std::fabs( (temperatureNumericalPtr[i] - temperatureTargetPtr[i])/temperatureTargetPtr[i] ) << std::endl << std::endl;
          }*/

          convergenceStatus = false;
        }
        //else
        //{
          // Don't modify guess if we are already converged
          //temperaturePtr[i] = temperatureLastGuessPtr[i];
        //}
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
