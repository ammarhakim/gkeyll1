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

    applyLeftEdge = applyRightEdge = false;
    std::string edgeStr = tbl.getString("edge");
    if (edgeStr == "lower")
      applyLeftEdge = true;
    else if (edgeStr == "upper")
      applyRightEdge = true;
    else if (edgeStr == "both")
      applyLeftEdge = applyRightEdge = true;

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

#ifdef HAVE_MPI
// barf if we this is called in parallel
    throw Lucee::Except("DistFuncReflectionBcUpdater does not work in parallel!");
#endif

    Lucee::Region<2, int> globalRgn = grid.getGlobalRegion();
    Lucee::FieldPtr<double> sknPtr = distf.createPtr(); // for skin-cell
    Lucee::FieldPtr<double> gstPtr = distf.createPtr(); // for ghost-cell

    Lucee::Matrix<double> nodeCoords(nodalBasis->getNumNodes(), 3);
    unsigned numNodes = nodalBasis->getNumNodes();

    if (applyRightEdge)
    {
      int ix = globalRgn.getUpper(0)-1; // right skin cell x index
      for (unsigned js=globalRgn.getUpper(1)-1, jg=0; js>=0; --js, ++jg)
      {
        nodalBasis->setIndex(ix, js);
        nodalBasis->getNodalCoordinates(nodeCoords);

        distf.setPtr(sknPtr, ix, js);
        distf.setPtr(gstPtr, ix+1, jg);

        for (unsigned k=0; k<numNodes; ++k)
        {
          if (std::fabs(nodeCoords(k,1)) < cutOffVel)
// copy data into ghost after rotating skin cell data by 180 degrees
            gstPtr[k] = sknPtr[rotMapRight[k]];
          else
// set to no inflow condition
            gstPtr[k] = 0.0;
        }
      }
    }
    
    if (applyLeftEdge)
    {
      int ix = globalRgn.getLower(0); // left skin cell x index
      for (unsigned js=globalRgn.getUpper(1)-1, jg=0; js>=0; --js, ++jg)
      {
        nodalBasis->setIndex(ix, js);
        nodalBasis->getNodalCoordinates(nodeCoords);

        distf.setPtr(sknPtr, ix, js);
        distf.setPtr(gstPtr, ix-1, jg);

        for (unsigned k=0; k<numNodes; ++k)
        {
          if (std::fabs(nodeCoords(k,1)) < cutOffVel)
// copy data into ghost after rotating skin cell data by 180 degrees
            gstPtr[k] = sknPtr[rotMapRight[k]];
          else
// set to no inflow condition
            gstPtr[k] = 0.0; 
        }
      }
    }

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
