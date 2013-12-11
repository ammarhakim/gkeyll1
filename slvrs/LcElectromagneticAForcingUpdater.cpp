/**
 * @file	LcElectromagneticAForcingUpdater.cpp
 *
 * @brief	Updater to compute A_parallel use inputs calculated in x-p space
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcElectromagneticAForcingUpdater.h>

namespace Lucee
{
// set id for module system
  const char *ElectromagneticAForcingUpdater::id = "ElectromagneticAForcingUpdater";

  ElectromagneticAForcingUpdater::ElectromagneticAForcingUpdater()
    : UpdaterIfc()
  {
  }

  void 
  ElectromagneticAForcingUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis");
    else
      throw Lucee::Except("ElectromagneticAForcingUpdater::readInput: Must specify element to use using 'basis'");

    if (tbl.hasNumber("kPerp"))
      kPerp = tbl.getNumber("kPerp");
    else
      throw Lucee::Except("ElectromagneticAForcingUpdater::readInput: Must specify kPerp");

    if (tbl.hasNumber("elcMass"))
      elcMass = tbl.getNumber("elcMass");
    else
      throw Lucee::Except("ElectromagneticAForcingUpdater::readInput: Must specify elcMass");

    if (tbl.hasNumber("ionMass"))
      ionMass = tbl.getNumber("ionMass");
    else
      throw Lucee::Except("ElectromagneticAForcingUpdater::readInput: Must specify ionMass");

    if (tbl.hasNumber("elcCharge"))
      elcCharge = tbl.getNumber("elcCharge");
    else
      throw Lucee::Except("ElectromagneticAForcingUpdater::readInput: Must specify elcCharge");

    if (tbl.hasNumber("ionCharge"))
      ionCharge = tbl.getNumber("ionCharge");
    else
      throw Lucee::Except("ElectromagneticAForcingUpdater::readInput: Must specify ionCharge");

    if (tbl.hasNumber("mu0"))
      mu0 = tbl.getNumber("mu0");
    else
      throw Lucee::Except("ElectromagneticAForcingUpdater::readInput: Must specify mu0");
  }

  void 
  ElectromagneticAForcingUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();
    // global region to update
    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion();

    Lucee::RowMajorSequencer<1> seq(globalRgn);
    seq.step(); // just to get to first index
    int idx[1];
    seq.fillWithIndex(idx);
    nodalBasis->setIndex(idx);
    
    int nlocal = nodalBasis->getNumNodes();
    int numQuadNodes = nodalBasis->getNumGaussNodes();

    // Get mass matrix and then copy to Eigen format
    Lucee::Matrix<double> massMatrixLucee(nlocal, nlocal);

    massMatrix = Eigen::MatrixXd(nlocal, nlocal);
    
    nodalBasis->getMassMatrix(massMatrixLucee);
    
    // Get interpolation matrix, gaussian quadrature points, and weights
    Lucee::Matrix<double> interpMatrixLucee(numQuadNodes, nlocal);

    Lucee::Matrix<double> gaussOrdinatesLucee(numQuadNodes, 3);

    std::vector<double> gaussWeights(numQuadNodes);

    // Allocate Eigen matrices
    Eigen::MatrixXd interpMatrix(numQuadNodes, nlocal);
    Eigen::MatrixXd gaussOrdinates(numQuadNodes, 3);

    // Get the interpolation matrix for the volume quadrature points
    nodalBasis->getGaussQuadData(interpMatrixLucee, gaussOrdinatesLucee, gaussWeights);

    copyLuceeToEigen(massMatrixLucee, massMatrix);
    copyLuceeToEigen(interpMatrixLucee, interpMatrix);
    copyLuceeToEigen(gaussOrdinatesLucee, gaussOrdinates);

    tripleProducts = std::vector<Eigen::MatrixXd>(nlocal);

    // Create and store the triple-product basis integrals
    for (int basisIndex = 0; basisIndex < nlocal; basisIndex++)
    {
      // Initialize matrices
      tripleProducts[basisIndex]  = Eigen::MatrixXd::Zero(nlocal, nlocal);
      for (int rowIndex = 0; rowIndex < nlocal; rowIndex++)
      {
        for (int colIndex = 0; colIndex < nlocal; colIndex++)
        {
          double integralResult = 0.0;

          for (int gaussNodeIndex = 0; gaussNodeIndex < numQuadNodes; gaussNodeIndex++)
          {
            integralResult += gaussWeights[gaussNodeIndex]*interpMatrix(gaussNodeIndex, basisIndex)*
              interpMatrix(gaussNodeIndex, rowIndex)*interpMatrix(gaussNodeIndex, colIndex);
          }

          tripleProducts[basisIndex](rowIndex, colIndex) = integralResult;
        }
      }
    }

    // Store inverse of mass matrix
    massMatrixInv = massMatrix.inverse();
  }

  Lucee::UpdaterStatus 
  ElectromagneticAForcingUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<1>& grid
      = this->getGrid<Lucee::StructuredGridBase<1> >();

    // Electron and ion densities
    const Lucee::Field<1, double>& mom0ElcIn = this->getInp<Lucee::Field<1, double> >(0);
    const Lucee::Field<1, double>& mom0IonIn = this->getInp<Lucee::Field<1, double> >(1);
    const Lucee::Field<1, double>& mom1ElcIn = this->getInp<Lucee::Field<1, double> >(2);
    const Lucee::Field<1, double>& mom1IonIn = this->getInp<Lucee::Field<1, double> >(3);
    const Lucee::Field<1, double>& aForcedIn = this->getInp<Lucee::Field<1, double> >(4);
    Lucee::Field<1, double>& aOut = this->getOut<Lucee::Field<1, double> >(0);

    int nlocal = nodalBasis->getNumNodes();

    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion();

    Lucee::ConstFieldPtr<double> mom0ElcPtr = mom0ElcIn.createConstPtr();
    Lucee::ConstFieldPtr<double> mom0IonPtr = mom0IonIn.createConstPtr();
    Lucee::ConstFieldPtr<double> mom1ElcPtr = mom1ElcIn.createConstPtr();
    Lucee::ConstFieldPtr<double> mom1IonPtr = mom1IonIn.createConstPtr();
    Lucee::ConstFieldPtr<double> aForcedPtr = aForcedIn.createConstPtr();
    Lucee::FieldPtr<double> aPtr = aOut.createPtr();

    aPtr = 0.0;
    
    // Loop over all cells
    for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
    {
      // Set inputs
      mom0ElcIn.setPtr(mom0ElcPtr, ix);
      mom0IonIn.setPtr(mom0IonPtr, ix);
      mom1ElcIn.setPtr(mom1ElcPtr, ix);
      mom1IonIn.setPtr(mom1IonPtr, ix);
      aForcedIn.setPtr(aForcedPtr, ix);
      // Set outputs
      aOut.setPtr(aPtr, ix);

      // Part 1 are terms linear in basis functions
      Eigen::VectorXd rhsVecPart1(nlocal);
      // Part 2 are terms quadratic in basis functions
      Eigen::VectorXd rhsVecPart2(nlocal);
      Eigen::VectorXd lhsVec(nlocal);
      // Vector of forced A_par components
      Eigen::VectorXd aForcedVec(nlocal);
      
      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
      {
        rhsVecPart1(componentIndex) = mu0/(kPerp*kPerp)*( 
            ionCharge*mom1IonPtr[componentIndex]/(ionMass*ionMass) + 
            elcCharge*mom1ElcPtr[componentIndex]/(elcMass*elcMass) );
        rhsVecPart2(componentIndex) = -mu0/(kPerp*kPerp)*(
            elcCharge*elcCharge*mom0ElcPtr[componentIndex]/(elcMass*elcMass) + 
            ionCharge*ionCharge*mom0IonPtr[componentIndex]/(ionMass*ionMass) );
        aForcedVec(componentIndex) = aForcedPtr[componentIndex];
      }

      Eigen::MatrixXd rhsPart2Projections(nlocal, nlocal);

      for (int basisIndex = 0; basisIndex < nlocal; basisIndex++)
      {
        Eigen::VectorXd rhsPart2ProjectionsSingle = tripleProducts[basisIndex]*rhsVecPart2;
        // Store components into rows of aProjections matrix
        for (int colIndex = 0; colIndex < nlocal; colIndex++)
          rhsPart2Projections(basisIndex, colIndex) = rhsPart2ProjectionsSingle(colIndex);
      }

      // Compute A(x)*n(x)*phi_n(x), where phi_n(x) is a basis function
      // Then combine with part 1 basis function integrals
      Eigen::VectorXd rhsComplete = massMatrix*rhsVecPart1 + rhsPart2Projections*aForcedVec;
      Eigen::VectorXd aWeights = massMatrixInv*rhsComplete;

      /*
      // Extra step to rescale A^{n+1} using A_out to test stability
      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
      {
        double aScaledValue = 1.0/(1.0 + mu0*elcCharge*elcCharge*mom0ElcPtr[componentIndex]/(elcMass*elcMass*kPerp*kPerp)
            + mu0*ionCharge*ionCharge*mom0IonPtr[componentIndex]/(ionMass*ionMass*kPerp*kPerp) )*aWeights(componentIndex);
        aWeights(componentIndex) = aScaledValue;
      }*/

      // Copy into output pointer
      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
        aPtr[componentIndex] = aWeights(componentIndex);
    }

    return Lucee::UpdaterStatus();
  }

  void
  ElectromagneticAForcingUpdater::declareTypes()
  {
    // takes mom0Elc, mom0Ion, mom1Elc, mom1Ion (in p-space)
    // last input: A_par (forced term)
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // returns one output: A(x)
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }

  void
  ElectromagneticAForcingUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}
