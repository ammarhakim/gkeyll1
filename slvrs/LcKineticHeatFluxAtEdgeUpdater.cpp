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
#include <LcKineticHeatFluxAtEdgeUpdater.h>
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

    if (tbl.hasNumber("electronMass"))
      elcMass = tbl.getNumber("electronMass");
    else
      throw Lucee::Except("KineticHeatFluxAtEdgeUpdater::readInput: Must specify electronMass");

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

    // Check to see if we need to compute the sheat power transmission coefficient
    computeSheathCoefficient = false;
    if (tbl.hasBool("computeSheathCoefficient"))
      computeSheathCoefficient = tbl.getBool("computeSheathCoefficient");
  }

  void 
  KineticHeatFluxAtEdgeUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus 
  KineticHeatFluxAtEdgeUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<1>& grid
      = this->getGrid<Lucee::StructuredGridBase<1> >();

    // potential
    const Lucee::Field<1, double>& phiIn = this->getInp<Lucee::Field<1, double> >(0);
    // mom1 and mom3 at left and right edges
    const Lucee::DynVector<double>& momentsAtEdgesElcIn = this->getInp<Lucee::DynVector<double> >(1);
    const Lucee::DynVector<double>& momentsAtEdgesIonIn = this->getInp<Lucee::DynVector<double> >(2);
    // Returns heat flux vs time as a dynvector
    Lucee::DynVector<double>& qVsTime = this->getOut<Lucee::DynVector<double> >(0);

    int nlocal = nodalBasis->getNumNodes();

    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion();

    Lucee::ConstFieldPtr<double> phiPtr = phiIn.createConstPtr();

    std::vector<double> momentsAtEdgesElc = momentsAtEdgesElcIn.getLastInsertedData();
    std::vector<double> momentsAtEdgesIon = momentsAtEdgesIonIn.getLastInsertedData();
    
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

    // Find value of the following input fields at the right-most edge of the domain
    phiIn.setPtr(phiPtr, globalRgn.getUpper(0)-1);
    
    double ionParallelHeatFluxRight = 0.5*ionMass*momentsAtEdgesIon[7] + 
      momentsAtEdgesIon[5]*ELEMENTARY_CHARGE*phiPtr[nlocal-1];
    double ionHeatFluxRight = ionParallelHeatFluxRight +
      momentsAtEdgesIon[5]*ELEMENTARY_CHARGE*tPerpIon;

    double elcParallelHeatFluxRight = 0.5*elcMass*momentsAtEdgesElc[7] - 
      momentsAtEdgesElc[5]*ELEMENTARY_CHARGE*phiPtr[nlocal-1];
    double elcHeatFluxRight = elcParallelHeatFluxRight + 
      momentsAtEdgesElc[5]*ELEMENTARY_CHARGE*tPerpElc;

    // Find value of the following input fields at the left-most edge of the domain
    phiIn.setPtr(phiPtr, globalRgn.getLower(0));

    double ionParallelHeatFluxLeft = 0.5*ionMass*momentsAtEdgesIon[3] + 
      momentsAtEdgesIon[1]*ELEMENTARY_CHARGE*phiPtr[0];
    double ionHeatFluxLeft = ionParallelHeatFluxLeft +
      momentsAtEdgesIon[1]*ELEMENTARY_CHARGE*tPerpIon;

    double elcParallelHeatFluxLeft = 0.5*elcMass*momentsAtEdgesElc[3] - 
      momentsAtEdgesElc[1]*ELEMENTARY_CHARGE*phiPtr[0];
    double elcHeatFluxLeft = elcParallelHeatFluxLeft + 
      momentsAtEdgesElc[1]*ELEMENTARY_CHARGE*tPerpElc;

    std::vector<double> data(6);
    data[0] = ionHeatFluxRight + elcHeatFluxRight;
    data[1] = ionHeatFluxRight;
    data[2] = elcHeatFluxRight;
    data[3] = -(ionHeatFluxLeft + elcHeatFluxLeft);
    data[4] = -ionHeatFluxLeft;
    data[5] = -elcHeatFluxLeft;
    
    qVsTime.appendData(t, data);

    if (computeSheathCoefficient == true)
    {
      Lucee::DynVector<double>& sheathCoefficientVsTime = this->getOut<Lucee::DynVector<double> >(1);
      
      std::vector<double> sheathData(1);

      double elcParallelTempRight = -elcMass*momentsAtEdgesElc[5]*momentsAtEdgesElc[5]/
        (momentsAtEdgesElc[4]*momentsAtEdgesElc[4]) + elcMass*momentsAtEdgesElc[6]/momentsAtEdgesElc[4];
            
      double kTe0 = 1.0/3.0*elcParallelTempRight + 2.0/3.0*ELEMENTARY_CHARGE*tPerpElc;
      sheathData[0] = data[0]/( kTe0*momentsAtEdgesElc[5] );
      
      /*
      double ionParallelTempRight = -ionMass*momentsAtEdgesIon[5]*momentsAtEdgesIon[5]/
        (momentsAtEdgesIon[4]*momentsAtEdgesIon[4]) + ionMass*momentsAtEdgesIon[6]/momentsAtEdgesIon[4];
      double elcParallelTempRight = -elcMass*momentsAtEdgesElc[5]*momentsAtEdgesElc[5]/
        (momentsAtEdgesElc[4]*momentsAtEdgesElc[4]) + elcMass*momentsAtEdgesElc[6]/momentsAtEdgesElc[4];


      sheathData[0] = ionParallelHeatFluxRight/(ionParallelTempRight*momentsAtEdgesIon[5]);
      sheathData[1] = elcParallelHeatFluxRight/(elcParallelTempRight*momentsAtEdgesElc[5]);

      double ionParallelTempLeft = -ionMass*momentsAtEdgesIon[1]*momentsAtEdgesIon[1]/
        (momentsAtEdgesIon[0]*momentsAtEdgesIon[0]) + ionMass*momentsAtEdgesIon[2]/momentsAtEdgesIon[0];
      double elcParallelTempLeft = -elcMass*momentsAtEdgesElc[1]*momentsAtEdgesElc[1]/
        (momentsAtEdgesElc[0]*momentsAtEdgesElc[0]) + elcMass*momentsAtEdgesElc[2]/momentsAtEdgesElc[0];

      sheathData[2] = ionParallelHeatFluxLeft/(ionParallelTempLeft*momentsAtEdgesIon[1]);
      sheathData[3] = elcParallelHeatFluxLeft/(elcParallelTempLeft*momentsAtEdgesElc[1]);*/

      sheathCoefficientVsTime.appendData(t, sheathData);
    }

    return Lucee::UpdaterStatus();
  }

  void
  KineticHeatFluxAtEdgeUpdater::declareTypes()
  {
    // takes phi(x) + dynvectors containing moments at edges of domain
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::DynVector<double>));
    this->appendInpVarType(typeid(Lucee::DynVector<double>));
    // returns one output: dynvector of heat flux at edge vs time
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
    // optional second output: dynvector of sheat power transmission
    // coefficient vs time
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
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
