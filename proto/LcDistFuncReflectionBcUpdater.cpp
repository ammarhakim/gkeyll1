/**
 * @file	LcDistFuncReflectionBcUpdater.cpp
 *
 * @brief	Project a function of a basis functions.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDistFuncReflectionBcUpdater.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <limits>
#include <vector>

namespace Lucee
{
  static const unsigned LC_LOWER_EDGE = 0;
  static const unsigned LC_UPPER_EDGE = 1;
  static const unsigned LC_BOTH_EDGES = 2;

  const char *DistFuncReflectionBcUpdater::id = "DistFuncReflectionBc";


  DistFuncReflectionBcUpdater::DistFuncReflectionBcUpdater()
  {
  }

  void
  DistFuncReflectionBcUpdater::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis");
    else
      throw Lucee::Except("DistFuncReflectionBcUpdater::readInput: Must specify element to use using 'basis'");

    edge = LC_BOTH_EDGES; // by default apply to both sides of domain
    std::string edgeStr = tbl.getString("edge");
    if (edgeStr == "lower")
      edge = LC_LOWER_EDGE;
    else if (edgeStr == "upper")
      edge = LC_UPPER_EDGE;

    cutOffVel = std::numeric_limits<double>::max(); // cut-off at infinity
    if (tbl.hasNumber("cutOffVelocity"))
      cutOffVel = tbl.getNumber("cutOffVelocity");
  }

  void
  DistFuncReflectionBcUpdater::initialize()
  {
    UpdaterIfc::initialize();

    unsigned numNodes = nodalBasis->getNumNodes();
    std::vector<unsigned> yRef(numNodes), xRef(numNodes);

// right side reflection mapping
    std::vector<unsigned> rotMapRight(numNodes);
    nodalBasis->getUpperReflectingBcMapping(0, yRef);
    nodalBasis->getLowerReflectingBcMapping(1, xRef);
    for (unsigned i=0; i<numNodes; ++i)
      rotMapRight[i] = xRef[yRef[i]];

// left side reflection mapping
    std::vector<unsigned> rotMapLeft(numNodes);
    nodalBasis->getLowerReflectingBcMapping(0, yRef);
    nodalBasis->getUpperReflectingBcMapping(1, xRef);
    for (unsigned i=0; i<numNodes; ++i)
      rotMapLeft[i] = xRef[yRef[i]];
  }

  Lucee::UpdaterStatus
  DistFuncReflectionBcUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();
    Lucee::Field<2, double>& distf = this->getOut<Lucee::Field<2, double> >(0);

    return Lucee::UpdaterStatus();
  }

  void
  DistFuncReflectionBcUpdater::declareTypes()
  {
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }

  void
  DistFuncReflectionBcUpdater::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
    Lucee::UpdaterIfc::appendLuaCallableMethods(lfm);

    lfm.appendFunc("setCutOffVelocity", luaSetCutOffVelocity);
  }

  int
  DistFuncReflectionBcUpdater::luaSetCutOffVelocity(lua_State *L)
  {
    DistFuncReflectionBcUpdater *up
      = Lucee::PointerHolder<DistFuncReflectionBcUpdater>::getObjAsDerived(L);
    double cv = lua_tonumber(L, 2);
    up->setCutOffVelocity(cv);
    return 0;
  }
}
