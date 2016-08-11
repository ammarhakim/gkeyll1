/**
 * @file	LcSOLMaxwellianParameterNewtonCalc.cpp
 *
 * @brief	Performs an iteration of Newton's method to calculate next parameters to use as
 * inputs for SOLMaxwellianAtNodeCalc
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLMaxwellianParameterNewtonCalc.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  const char *SOLMaxwellianParameterNewtonCalc::id = "SOLMaxwellianParameterNewtonCalc";

  SOLMaxwellianParameterNewtonCalc::SOLMaxwellianParameterNewtonCalc()
  {
  }

  void
  SOLMaxwellianParameterNewtonCalc::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("SOLMaxwellianParameterNewtonCalc::readInput: Must specify element to use using 'basis3d'");

    if (tbl.hasNumber("epsilonU"))
      epsilonU = tbl.getNumber("epsilonU");
    else 
      throw Lucee::Except("SOLMaxwellianParameterNewtonCalc::readInput: Must specify epsilonU to use using 'epsilonU'");

    if (tbl.hasNumber("epsilonT"))
      epsilonT = tbl.getNumber("epsilonT");
    else 
      throw Lucee::Except("SOLMaxwellianParameterNewtonCalc::readInput: Must specify epsilonT to use using 'epsilonT'");
  }

  void
  SOLMaxwellianParameterNewtonCalc::initialize()
  {
    UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus
  SOLMaxwellianParameterNewtonCalc::update(double t)
  {
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    const Lucee::Field<3, double>& mom1dir3TargetIn = this->getInp<Lucee::Field<3, double> >(0);
    const Lucee::Field<3, double>& temperatureTargetIn = this->getInp<Lucee::Field<3, double> >(1);
    const Lucee::Field<3, double>& mom1dir3NumericalIn = this->getInp<Lucee::Field<3, double> >(2);
    const Lucee::Field<3, double>& temperatureNumericalIn = this->getInp<Lucee::Field<3, double> >(3);
    const Lucee::Field<3, double>& mom1dir3IncrementOneIn = this->getInp<Lucee::Field<3, double> >(4);
    const Lucee::Field<3, double>& temperatureIncrementOneIn = this->getInp<Lucee::Field<3, double> >(5);
    const Lucee::Field<3, double>& mom1dir3IncrementTwoIn = this->getInp<Lucee::Field<3, double> >(6);
    const Lucee::Field<3, double>& temperatureIncrementTwoIn = this->getInp<Lucee::Field<3, double> >(7);
    const Lucee::Field<3, double>& mom1dir3InputIn = this->getInp<Lucee::Field<3, double> >(8);
    const Lucee::Field<3, double>& temperatureInputIn = this->getInp<Lucee::Field<3, double> >(9);
    // Output: Next guess increments for SOLMaxwellianAtNodeCalc
    Lucee::Field<3, double>& mom1dir3Out = this->getOut<Lucee::Field<3, double> >(0);
    Lucee::Field<3, double>& temperatureOut = this->getOut<Lucee::Field<3, double> >(1);
    Lucee::Field<3, double>& slopeOut = this->getOut<Lucee::Field<3, double> >(2);

    Lucee::ConstFieldPtr<double> mom1dir3TargetPtr = mom1dir3TargetIn.createConstPtr();
    Lucee::ConstFieldPtr<double> temperatureTargetPtr = temperatureTargetIn.createConstPtr();

    Lucee::ConstFieldPtr<double> mom1dir3NumericalPtr = mom1dir3NumericalIn.createConstPtr();
    Lucee::ConstFieldPtr<double> temperatureNumericalPtr = temperatureNumericalIn.createConstPtr();

    // Numerical values for epsilon increment in mom1dir3
    Lucee::ConstFieldPtr<double> mom1dir3IncrementOnePtr = mom1dir3IncrementOneIn.createConstPtr();
    Lucee::ConstFieldPtr<double> temperatureIncrementOnePtr = temperatureIncrementOneIn.createConstPtr();

    // Numerical values for epsilon increment in temperature
    Lucee::ConstFieldPtr<double> mom1dir3IncrementTwoPtr = mom1dir3IncrementTwoIn.createConstPtr();
    Lucee::ConstFieldPtr<double> temperatureIncrementTwoPtr = temperatureIncrementTwoIn.createConstPtr();

    // Input parameters to compute Maxwellian at this step
    Lucee::ConstFieldPtr<double> mom1dir3InputPtr = mom1dir3InputIn.createConstPtr();
    Lucee::ConstFieldPtr<double> temperatureInputPtr = temperatureInputIn.createConstPtr();
    
    // Output to be used in next maxwellian calculation
    Lucee::FieldPtr<double> mom1dir3Ptr = mom1dir3Out.createPtr();
    Lucee::FieldPtr<double> temperaturePtr = temperatureOut.createPtr();
    Lucee::FieldPtr<double> slopePtr = slopeOut.createPtr();

    unsigned nlocal3d = nodalBasis3d->getNumNodes();

    int idx[3];
    Lucee::Region<3, int> localRgn = grid.getLocalRegion();
    Lucee::RowMajorSequencer<3> seq(localRgn);

    const double functionTol = 1e-4;
    const double rootTol = 1e-8;
    bool rootConvergenceStatus = true;
    bool functionConvergenceStatus = true;

    Eigen::Matrix2d jacobianMatrix;
    Eigen::Vector2d fVec;
    Eigen::Vector2d solutionVec;
    Eigen::Vector2d gradVec;

    // Loop over local region
    while (seq.step())
    {
      seq.fillWithIndex(idx);

      // Set input fields to right location
      mom1dir3TargetIn.setPtr(mom1dir3TargetPtr, idx);
      mom1dir3NumericalIn.setPtr(mom1dir3NumericalPtr, idx);
      mom1dir3IncrementOneIn.setPtr(mom1dir3IncrementOnePtr, idx);
      mom1dir3IncrementTwoIn.setPtr(mom1dir3IncrementTwoPtr, idx);
      mom1dir3InputIn.setPtr(mom1dir3InputPtr, idx);
      
      temperatureTargetIn.setPtr(temperatureTargetPtr, idx);
      temperatureNumericalIn.setPtr(temperatureNumericalPtr, idx);
      temperatureIncrementOneIn.setPtr(temperatureIncrementOnePtr, idx);
      temperatureIncrementTwoIn.setPtr(temperatureIncrementTwoPtr, idx);
      temperatureInputIn.setPtr(temperatureInputPtr, idx);
      // Set output fields to right location
      mom1dir3Out.setPtr(mom1dir3Ptr, idx);
      temperatureOut.setPtr(temperaturePtr, idx);
      slopeOut.setPtr(slopePtr, idx);

      for (int i = 0; i < nlocal3d; i++)
      {
        // Compute jacobian matrix using forward differences
        jacobianMatrix(0,0) = (mom1dir3IncrementOnePtr[i]-mom1dir3NumericalPtr[i])/( epsilonU*std::max(fabs(mom1dir3InputPtr[i]), 1.0) );
        jacobianMatrix(0,1) = 0.0;//(mom1dir3IncrementTwoPtr[i]-mom1dir3NumericalPtr[i])/(epsilonT*fabs(temperatureInputPtr[i]));
        jacobianMatrix(1,0) = (temperatureIncrementOnePtr[i]-temperatureNumericalPtr[i])/( epsilonU*std::max(fabs(mom1dir3InputPtr[i]), 1.0) );
        jacobianMatrix(1,1) = (temperatureIncrementTwoPtr[i]-temperatureNumericalPtr[i])/(epsilonT*fabs(temperatureInputPtr[i]));

        // Fill out fVec with a minus sign
        fVec(0) = -(mom1dir3NumericalPtr[i] - mom1dir3TargetPtr[i]);
        fVec(1) = -(temperatureNumericalPtr[i] - temperatureTargetPtr[i]);

        // Solve for increment
        solutionVec = jacobianMatrix.colPivHouseholderQr().solve(fVec);

        // Set next guess values
        mom1dir3Ptr[i] = solutionVec(0);
        // Limit next iteration guess on mom1dir3Ptr[i] to 5x the magnitude of mom1dir3Ptr
        //if (std::fabs(solutionVec(0)) > 5*std::fabs(mom1dir3InputPtr[i]))
        //  mom1dir3Ptr[i] = mom1dir3InputPtr[i] + 5*solutionVec(0)/std::fabs(solutionVec(0))*mom1dir3InputPtr[i];

        temperaturePtr[i] = solutionVec(1);
        // If next temperature guess is to be negative, just take half the value
        //if (temperaturePtr[i] <= 0.0)
        //  temperaturePtr[i] = temperatureInputPtr[i]*0.5;
        
        // Compute gradient of 0.5*fVec.dot(fVec)
        gradVec = jacobianMatrix.transpose()*fVec;
        // Store "slope" variable gradfVec dot p
        slopePtr[i] = gradVec.dot(solutionVec);

        // Check convergence in increment size
        if (std::fabs(solutionVec(0)) > std::fabs(mom1dir3InputPtr[i]*rootTol))
          rootConvergenceStatus = false;

        if (std::fabs(solutionVec(1)) > std::fabs(temperatureInputPtr[i]*rootTol))
          rootConvergenceStatus = false;

        // Check convergence in function
        if ( std::fabs(mom1dir3NumericalPtr[i] - mom1dir3TargetPtr[i]) > std::fabs(mom1dir3TargetPtr[i]*functionTol)
             || std::fabs(temperatureNumericalPtr[i] - temperatureTargetPtr[i]) > std::fabs(temperatureTargetPtr[i]*functionTol) )
        {
          std::cout << "idx[" << idx[0] << "," << idx[1] << "," << idx[2] << "]" << std::endl;
          std::cout << "mom1dir3TargetPtr " << i << " = " << mom1dir3TargetPtr[i] << std::endl;
          std::cout << "mom1dir3NumericalPtr " << i << " = " << mom1dir3NumericalPtr[i] << std::endl;
          std::cout << "mom1dir3IncrementOnePtr " << i << " = " << mom1dir3IncrementOnePtr[i] << std::endl;
          std::cout << "mom1dir3IncrementTwoPtr " << i << " = " << mom1dir3IncrementTwoPtr[i] << std::endl;
          std::cout << "mom1dir3LastGuessPtr " << i << " = " << mom1dir3InputPtr[i] << std::endl;
          std::cout << "mom1dirNextGuess " << i << " = " << mom1dir3Ptr[i] << std::endl;
          std::cout << "temperatureTargetPtr " << i << " = " << temperatureTargetPtr[i] << std::endl;
          std::cout << "temperatureNumericalPtr " << i << " = " << temperatureNumericalPtr[i] << std::endl;
          std::cout << "temperatureIncrementOnePtr " << i << " = " << temperatureIncrementOnePtr[i] << std::endl;
          std::cout << "temperatureIncrementTwoPtr " << i << " = " << temperatureIncrementTwoPtr[i] << std::endl;
          std::cout << "temperatureLastGuessPtr " << i << " = " << temperatureInputPtr[i] << std::endl;
          std::cout << "temperatureNextGuess " << i << " = " << temperaturePtr[i] << std::endl;
          std::cout << "relErrorU = " << std::fabs( (mom1dir3NumericalPtr[i] - mom1dir3TargetPtr[i])/mom1dir3TargetPtr[i] ) << std::endl;
          std::cout << "relErrorT = " << std::fabs( (temperatureNumericalPtr[i] - temperatureTargetPtr[i])/temperatureTargetPtr[i] ) << std::endl;
          std::cout << "jacobianMatrix" << std::endl << jacobianMatrix << std::endl;
          std::cout << "solutionVec" << std::endl << solutionVec << std::endl << std::endl;
          functionConvergenceStatus = false;
        }
      }
    }

    // Insert a mpi call to "and" all the processes

    // Returns true if function or root has converged
    return Lucee::UpdaterStatus((functionConvergenceStatus || rootConvergenceStatus), 0.0);
  }

  void
  SOLMaxwellianParameterNewtonCalc::declareTypes()
  {
    // desired <v> = u*n
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // desired temperature (in joules)
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // numerical <v> = u*n
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // numerical temperature (in joules)
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // numerical <v> = u*n for epsilon increment in mom1dir3
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // numerical temperature (in joules) for epsilon increment in mom1dir3
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // numerical <v> = u*n for increment in temperature
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // numerical temperature (in joules) for increment in temperature
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // input <v> for this step to compute numerical maxwellian
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // input mom1dir3 for this step to compute numerical mawellian
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // next guess increment <v> = u*n
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
    // next guess increment temperature (in joules)
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
    // slope variable
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
  }
}
