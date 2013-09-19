/**
 * @file	LcHeatFluxAtEdgeUpdater.cpp
 *
 * @brief	Updater to compute heat flux at right-most edge in domain
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAlignedRectCoordSys.h>
#include <LcField.h>
#include <LcLinAlgebra.h>
#include <LcMathLib.h>
#include <LcHeatFluxAtEdgeUpdater.h>
#include <LcStructuredGridBase.h>
#include <LcMathPhysConstants.h>
// For function handle:
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
// set id for module system
  const char *HeatFluxAtEdgeUpdater::id = "HeatFluxAtEdgeUpdater";

  HeatFluxAtEdgeUpdater::HeatFluxAtEdgeUpdater()
    : UpdaterIfc()
  {
  }

  void 
  HeatFluxAtEdgeUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis");
    else
      throw Lucee::Except("HeatFluxAtEdgeUpdater::readInput: Must specify element to use using 'basis1d'");

    if (tbl.hasNumber("ionMass"))
      ionMass = tbl.getNumber("ionMass");
    else
      throw Lucee::Except("HeatFluxAtEdgeUpdater::readInput: Must specify ionMass");

    fnProvided = false;
    if (tbl.hasFunction("tPerpProfile"))
    {
      fnProvided = true;
      fnRef = tbl.getFunctionRef("tPerpProfile");
    }
    else
    {
      if (tbl.hasNumber("tPerp"))
        tPerp = tbl.getNumber("tPerp");
      else
        throw Lucee::Except("HeatFluxAtEdgeUpdater::readInput: Must specify tPerp");
    }
  }

  void 
  HeatFluxAtEdgeUpdater::initialize()
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

    nodalBasis->getMassMatrix(massMatrixLucee);
    
    // Get interpolation matrix, gaussian quadrature points, and weights
    Lucee::Matrix<double> interpMatrixLucee(numQuadNodes, nlocal);

    Lucee::Matrix<double> gaussOrdinatesLucee(numQuadNodes, 3);

    gaussWeights = std::vector<double>(numQuadNodes);

    // Allocate Eigen matrices
    interpMatrix = Eigen::MatrixXd(numQuadNodes, nlocal);

    // Get the interpolation matrix for the volume quadrature points
    nodalBasis->getGaussQuadData(interpMatrixLucee, gaussOrdinatesLucee, gaussWeights);
    copyLuceeToEigen(interpMatrixLucee, interpMatrix);
  }

  Lucee::UpdaterStatus 
  HeatFluxAtEdgeUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<1>& grid
      = this->getGrid<Lucee::StructuredGridBase<1> >();

    // First three velocity moments + driftU
    const Lucee::Field<1, double>& vtSqIn = this->getInp<Lucee::Field<1, double> >(0);
    const Lucee::Field<1, double>& phiIn = this->getInp<Lucee::Field<1, double> >(1);
    // mom1 and mom3 at left and right edges
    const Lucee::DynVector<double>& momentsAtEdgesIn = this->getInp<Lucee::DynVector<double> >(2);
    // Returns heat flux vs time as a dynvector
    Lucee::DynVector<double>& qVsTime = this->getOut<Lucee::DynVector<double> >(0);

    int nlocal = nodalBasis->getNumNodes();

    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion();

    Lucee::ConstFieldPtr<double> vtSqPtr = vtSqIn.createConstPtr();
    Lucee::ConstFieldPtr<double> phiPtr = phiIn.createConstPtr();

    std::vector<double> momentsAtEdges = momentsAtEdgesIn.getLastInsertedData();
    
    // Find mean vtSq in entire domain
    double meanVtSq = 0.0;

    for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
    {
      // Set inputs
      vtSqIn.setPtr(vtSqPtr, ix);
      Eigen::VectorXd vtSqVec(nlocal);
      
      // Figure out vtSq(x) at quadrature points in the cell
      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
        vtSqVec(componentIndex) = vtSqPtr[componentIndex];
      Eigen::VectorXd vtSqAtQuadPoints = interpMatrix*vtSqVec;

      for (int componentIndex = 0; componentIndex < vtSqAtQuadPoints.rows(); componentIndex++)
        meanVtSq += gaussWeights[componentIndex]*vtSqAtQuadPoints(componentIndex);
    }

    // Divide by length of domain
    // Consider using grid.getNumCells(0)
    meanVtSq = meanVtSq/(grid.getDx(0)*(globalRgn.getUpper(0)-globalRgn.getLower(0)));

    // Find value of the following input fields at the very last edge of the domain
    phiIn.setPtr(phiPtr, globalRgn.getUpper(0)-1);

    double tPerpIon;
    double tPerpElc;

    // Figure out value of perpendicular electron and ion temperatures
    if (fnProvided == true)
    {
      std::vector<double> res(2);
      Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;

      evaluateFunction(*L, t, res);

      tPerpIon = res[0];
      tPerpElc = res[1];
    }
    else
    {
      tPerpIon = tPerp;
      tPerpElc = tPerp;
    }
    
    double ionHeatFlux = 0.5*ionMass*momentsAtEdges[3] + momentsAtEdges[2]*ELEMENTARY_CHARGE*(tPerpIon + phiPtr[nlocal-1]);
    // Te = average of Ti(x)
    double electronHeatFlux = (ionMass*meanVtSq + ELEMENTARY_CHARGE*tPerpElc)*momentsAtEdges[2];
    
    std::vector<double> data(3);
    data[0] = ionHeatFlux + electronHeatFlux;
    data[1] = ionHeatFlux;
    data[2] = electronHeatFlux;
    qVsTime.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  void
  HeatFluxAtEdgeUpdater::declareTypes()
  {
    // takes inputs vtSq + phi(x)
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // moments evaluated at edges input
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
    // returns one output: dynvector of heat flux at edge vs time
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }
  
  void
  HeatFluxAtEdgeUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }

  void
  HeatFluxAtEdgeUpdater::evaluateFunction(Lucee::LuaState& L, double tm,
    std::vector<double>& res)
  {
    // push function object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, fnRef);
    // push variables on stack
    lua_pushnumber(L, tm);
    // call function
    if (lua_pcall(L, 1, res.size(), 0) != 0)
    {
      Lucee::Except lce("HeatFluxAtEdgeUpdater::evaluateFunction: ");
      lce << "Problem evaluating function supplied as 'tPerpProfile' "
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
        throw Lucee::Except("HeatFluxAtEdgeUpdater::evaluateFunction: Return value not a number");
      res[res.size()+i] = lua_tonumber(L, i);
    }
    lua_pop(L, 1);
  }
}
