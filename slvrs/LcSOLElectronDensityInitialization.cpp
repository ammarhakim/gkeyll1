/**
 * @file	LcSOLElectronDensityInitialization.cpp
 *
 * @brief	Updater to compute phi using a fixed value of k_perp*rho_s
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLElectronDensityInitialization.h>
#include <LcMathPhysConstants.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <vector>

namespace Lucee
{
// set id for module system
  const char *SOLElectronDensityInitialization::id = "SOLElectronDensityInitialization";

  SOLElectronDensityInitialization::SOLElectronDensityInitialization()
    : UpdaterIfc()
  {
  }

  void 
  SOLElectronDensityInitialization::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);
    // get function to evaluate
    fnRef = tbl.getFunctionRef("evaluate");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis");
    else
      throw Lucee::Except("SOLElectronDensityInitialization::readInput: Must specify element to use using 'basis'");

    if (tbl.hasNumber("kPerpTimesRho"))
      kPerpTimesRho = tbl.getNumber("kPerpTimesRho");
    else
      throw Lucee::Except("SOLElectronDensityInitialization::readInput: Must specify kPerpTimesRho");

    if (tbl.hasNumber("Te0"))
      Te0 = tbl.getNumber("Te0");
    else
      throw Lucee::Except("SOLElectronDensityInitialization::readInput: Must specify Te0");
  }

  void 
  SOLElectronDensityInitialization::initialize()
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

    // Get interpolation matrix, gaussian quadrature points, and weights
    Lucee::Matrix<double> interpMatrixLucee(numQuadNodes, nlocal);
    Lucee::Matrix<double> gaussOrdinatesLucee(numQuadNodes, 3);
    gaussWeights = std::vector<double>(numQuadNodes);
    // Get the interpolation matrix for the volume quadrature points
    nodalBasis->getGaussQuadData(interpMatrixLucee, gaussOrdinatesLucee, gaussWeights);

    // Allocate Eigen matrices
    interpMatrix = Eigen::MatrixXd(numQuadNodes, nlocal);

    copyLuceeToEigen(interpMatrixLucee, interpMatrix);
  }

  Lucee::UpdaterStatus 
  SOLElectronDensityInitialization::update(double t)
  {
    const Lucee::StructuredGridBase<1>& grid
      = this->getGrid<Lucee::StructuredGridBase<1> >();

    // Electron and ion densities
    const Lucee::Field<1, double>& nIonIn = this->getInp<Lucee::Field<1, double> >(0);
    Lucee::Field<1, double>& nElcOut = this->getOut<Lucee::Field<1, double> >(0);

    int nlocal = nodalBasis->getNumNodes();

    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion();

    Lucee::ConstFieldPtr<double> nIonInPtr = nIonIn.createConstPtr();
    Lucee::FieldPtr<double> nElcOutPtr = nElcOut.createPtr();

    nElcOutPtr = 0.0;

    // get hold of Lua state object
    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;
    // to store function evaluation result
    std::vector<double> resultVector(1);

    int numCells = globalRgn.getUpper(0)-globalRgn.getLower(0);
    Eigen::VectorXd zPrev(numCells*nlocal);

    // Calculate integrated ion density
    double intIonDensity = 0.0;
    for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
    {
      nIonIn.setPtr(nIonInPtr, ix);
      Eigen::VectorXd nIonVec(nlocal);

      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
        nIonVec(componentIndex) = nIonInPtr[componentIndex];
      // Calculate nIon at quadrature points
      Eigen::VectorXd nIonAtQuadPoints = interpMatrix*nIonVec;

      for (int componentIndex = 0; componentIndex < nIonAtQuadPoints.rows(); componentIndex++)
        intIonDensity += gaussWeights[componentIndex]*nIonAtQuadPoints(componentIndex);
    }
    
    int maxIter = 10000;
    // To get the correct position to evaluate temperature function as a position of z
    Lucee::Matrix<double> nodeCoordsLucee(nlocal, 3);

    int iter = 0;
    double tol = 1e-14;
    double relError;

    // Solve for z by iteration (store temporary results in output structure)
    do {
      // Calculate <exp(z)>
      double intExpZ = 0.0;
      for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
      {
        nElcOut.setPtr(nElcOutPtr, ix);
        Eigen::VectorXd expZVec(nlocal);

        for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
          expZVec(componentIndex) = std::exp(nElcOutPtr[componentIndex]);
       
        Eigen::VectorXd expZAtQuadPoints = interpMatrix*expZVec;
        
        for (int componentIndex = 0; componentIndex < expZAtQuadPoints.rows(); componentIndex++)
          intExpZ += gaussWeights[componentIndex]*expZAtQuadPoints(componentIndex);
      }

      // For calculating <expZ>
      double intExpZAgain = 0.0;
      // Loop over all cells
      for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
      {
        // Set inputs
        nIonIn.setPtr(nIonInPtr, ix);
        // Set outputs
        nElcOut.setPtr(nElcOutPtr, ix);

        // Set nodal basis index
        nodalBasis->setIndex(ix);
        // Get nodal coordinates into Lucee matrix
        nodalBasis->getNodalCoordinates(nodeCoordsLucee);

        Eigen::VectorXd expZVec(nlocal);

        for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
        {
          // Evaluate supplied temperature function
          evaluateFunction(*L, t, nodeCoordsLucee(nodeIndex, 0), resultVector);
          double elcTemp = resultVector[0];
          double logArg = nIonInPtr[nodeIndex]/intIonDensity*
             (1-kPerpTimesRho*kPerpTimesRho*nElcOutPtr[nodeIndex]*elcTemp/Te0)*
             intExpZ;
          
          nElcOutPtr[nodeIndex] = std::log(logArg);

          // Copy data into vector to perform integral of <expZ>
          expZVec(nodeIndex) = logArg;
        }

        Eigen::VectorXd expZAtQuadPoints = interpMatrix*expZVec;
        for (int componentIndex = 0; componentIndex < expZAtQuadPoints.rows(); componentIndex++)
          intExpZAgain += gaussWeights[componentIndex]*expZAtQuadPoints(componentIndex);
      }

      // Compute <(n_i - n_e)*Te0/kPerpRho^2>
      double convergenceFactor = 0.0;
      for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
      {
        // Set inputs
        nIonIn.setPtr(nIonInPtr, ix);
        // Set outputs
        nElcOut.setPtr(nElcOutPtr, ix);

        Eigen::VectorXd convergenceVec(nlocal);

        for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
          convergenceVec(nodeIndex) = (nIonInPtr[nodeIndex]-intIonDensity*std::exp(nElcOutPtr[nodeIndex])/
              intExpZAgain)*Te0/(kPerpTimesRho*kPerpTimesRho);
        
        Eigen::VectorXd convergenceVecAtQuadPoints = interpMatrix*convergenceVec;
        for (int quadIndex = 0; quadIndex < convergenceVecAtQuadPoints.rows(); quadIndex++)
          convergenceFactor += gaussWeights[quadIndex]*convergenceVecAtQuadPoints(quadIndex);
      }

      // For calculating relative errors
      double maxRes = 0.0;
      double maxZ = 0.0;
      // Go through all output locations and add correcting factor
      for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
      {
        // Set outputs
        nElcOut.setPtr(nElcOutPtr, ix);
        // Set nodal basis index
        nodalBasis->setIndex(ix);
        // Get nodal coordinates into Lucee matrix
        nodalBasis->getNodalCoordinates(nodeCoordsLucee);

        for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
        {
          // Evaluate supplied temperature function
          evaluateFunction(*L, t, nodeCoordsLucee(nodeIndex, 0), resultVector);
          double elcTemp = resultVector[0];
          nElcOutPtr[nodeIndex] = nElcOutPtr[nodeIndex] - convergenceFactor/(intIonDensity*elcTemp);
          
          if (fabs(nElcOutPtr[nodeIndex]) > maxZ)
            maxZ = fabs(nElcOutPtr[nodeIndex]);
          // Calculate and possibly store error
          double residual = fabs(nElcOutPtr[nodeIndex] - zPrev((ix-globalRgn.getLower(0))*nlocal + nodeIndex));
          if (residual > maxRes)
            maxRes = residual;
          // Store results for comparison
          zPrev((ix-globalRgn.getLower(0))*nlocal + nodeIndex) = nElcOutPtr[nodeIndex];
        }
      }

      relError = maxRes/maxZ;
      std::cout << "iter " << iter << " relerr = " << relError << std::endl;
      iter++;
    } while (iter < maxIter && relError > tol);

    // Calculate <exp(z)> again
    double intExpZ = 0.0;
    for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
    {
      nElcOut.setPtr(nElcOutPtr, ix);
      Eigen::VectorXd expZVec(nlocal);

      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
        expZVec(componentIndex) = std::exp(nElcOutPtr[componentIndex]);
     
      Eigen::VectorXd expZAtQuadPoints = interpMatrix*expZVec;
      
      for (int componentIndex = 0; componentIndex < expZAtQuadPoints.rows(); componentIndex++)
        intExpZ += gaussWeights[componentIndex]*expZAtQuadPoints(componentIndex);
    }

    // Final loop: take z and compute n_e = <n_i>*exp(z)/<exp(z)>
    for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
    {
      nElcOut.setPtr(nElcOutPtr, ix);

      for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
        nElcOutPtr[nodeIndex] = intIonDensity*std::exp(nElcOutPtr[nodeIndex])/intExpZ;
    }

    return Lucee::UpdaterStatus();
  }

  void
  SOLElectronDensityInitialization::declareTypes()
  {
    // takes one input: n_i(x), a 1d DG field
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // returns one output: n_e(x), a 1d DG field
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }

  void
  SOLElectronDensityInitialization::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }

  void
  SOLElectronDensityInitialization::evaluateFunction(Lucee::LuaState& L, double tm,
    double positionValue, std::vector<double>& res)
  {
    // push function object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, fnRef);
    // push variables on stack
    lua_pushnumber(L, positionValue);
    lua_pushnumber(L, tm);
    // call function
    if (lua_pcall(L, 2, res.size(), 0) != 0)
    {
      Lucee::Except lce("SOLElectronDensityInitialization::evaluateFunction: ");
      lce << "Problem evaluating function supplied as 'evaluate' "
          << std::endl;
      std::string err(lua_tostring(L, -1));
      lua_pop(L, 1);
      lce << "[" << err << "]";
      throw lce;
    }
    // fetch results
    for (int i=-res.size(); i<0; ++i)
    {
      if (!lua_isnumber(L, i))
        throw Lucee::Except("SOLElectronDensityInitialization::evaluateFunction: Return value not a number");
      res[res.size()+i] = lua_tonumber(L, i);
    }
    lua_pop(L, 1);
  }
}
