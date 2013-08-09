/**
 * @file	LcKineticKineticHeatFluxAtEdgeUpdater.cpp
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
#include <LcKineticHeatFluxAtEdgeUpdater.h>
#include <LcStructuredGridBase.h>
#include <LcMathPhysConstants.h>
// For function handle:
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
// set id for module system
  const char *KineticHeatFluxAtEdgeUpdater::id = "KineticHeatFluxAtEdgeUpdater";

  KineticHeatFluxAtEdgeUpdater::KineticHeatFluxAtEdgeUpdater()
    : UpdaterIfc()
  {
  }

  void 
  KineticHeatFluxAtEdgeUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis");
    else
      throw Lucee::Except("KineticHeatFluxAtEdgeUpdater::readInput: Must specify element to use using 'basis1d'");

    if (tbl.hasNumber("ionMass"))
      ionMass = tbl.getNumber("ionMass");
    else
      throw Lucee::Except("KineticHeatFluxAtEdgeUpdater::readInput: Must specify ionMass");

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
        throw Lucee::Except("KineticHeatFluxAtEdgeUpdater::readInput: Must specify tPerp");
    }
  }

  void 
  KineticHeatFluxAtEdgeUpdater::initialize()
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
  KineticHeatFluxAtEdgeUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<1>& grid
      = this->getGrid<Lucee::StructuredGridBase<1> >();

    // First three velocity moments + driftU
    const Lucee::Field<1, double>& mom1ElcIn = this->getInp<Lucee::Field<1, double> >(0);
    const Lucee::Field<1, double>& mom1IonIn = this->getInp<Lucee::Field<1, double> >(1);
    const Lucee::Field<1, double>& mom3ElcIn = this->getInp<Lucee::Field<1, double> >(2);
    const Lucee::Field<1, double>& mom3IonIn = this->getInp<Lucee::Field<1, double> >(3);
    const Lucee::Field<1, double>& phiIn     = this->getInp<Lucee::Field<1, double> >(4);
    // Returns heat flux vs time as a dynvector
    Lucee::DynVector<double>& qVsTime = this->getOut<Lucee::DynVector<double> >(0);

    int nlocal = nodalBasis->getNumNodes();

    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion();

    Lucee::ConstFieldPtr<double> mom1ElcPtr = mom1ElcIn.createConstPtr();
    Lucee::ConstFieldPtr<double> mom1IonPtr = mom1IonIn.createConstPtr();
    Lucee::ConstFieldPtr<double> mom3ElcPtr = mom3ElcIn.createConstPtr();
    Lucee::ConstFieldPtr<double> mom3IonPtr = mom3IonIn.createConstPtr();
    Lucee::ConstFieldPtr<double> phiPtr = phiIn.createConstPtr();

    // Find value of the following input fields at the very last edge of the domain
    mom1ElcIn.setPtr(mom1ElcPtr, globalRgn.getUpper(0)-1);
    mom1IonIn.setPtr(mom1IonPtr, globalRgn.getUpper(0)-1);
    mom3ElcIn.setPtr(mom3ElcPtr, globalRgn.getUpper(0)-1);
    mom3IonIn.setPtr(mom3IonPtr, globalRgn.getUpper(0)-1);
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
    
    double ionHeatFlux = 0.5*ionMass*mom3IonPtr[nlocal-1] + mom1IonPtr[nlocal-1]*ELEMENTARY_CHARGE*(tPerpIon + phiPtr[nlocal-1]);
    double electronHeatFlux = 0.5*ELECTRON_MASS*mom3ElcPtr[nlocal-1] + mom1ElcPtr[nlocal-1]*ELEMENTARY_CHARGE*(tPerpElc + phiPtr[nlocal-1]);
    
    std::vector<double> data(3);
    data[0] = ionHeatFlux + electronHeatFlux;
    data[1] = ionHeatFlux;
    data[2] = electronHeatFlux;
    qVsTime.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  void
  KineticHeatFluxAtEdgeUpdater::declareTypes()
  {
    // takes five inputs (<v>,<v^3>) elc moments + same for ions + phi(x)
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // returns one output: dynvector of heat flux at edge vs time
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }
  
  void
  KineticHeatFluxAtEdgeUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }

  void
  KineticHeatFluxAtEdgeUpdater::evaluateFunction(Lucee::LuaState& L, double tm,
    std::vector<double>& res)
  {
    // push function object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, fnRef);
    // push variables on stack
    lua_pushnumber(L, tm);
    // call function
    if (lua_pcall(L, 1, res.size(), 0) != 0)
    {
      Lucee::Except lce("KineticHeatFluxAtEdgeUpdater::evaluateFunction: ");
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
        throw Lucee::Except("KineticHeatFluxAtEdgeUpdater::evaluateFunction: Return value not a number");
      res[res.size()+i] = lua_tonumber(L, i);
    }
    lua_pop(L, 1);
  }
}
