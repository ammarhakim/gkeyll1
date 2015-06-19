/**
 * @file	LcDistFuncReflectionBcUpdater.cpp
 *
 * @brief	Applies particle refection BCs to (electron) distribution function.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcReflectionBoundaryCondition.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  static const unsigned LC_LOWER_EDGE = 0;
  static const unsigned LC_UPPER_EDGE = 1;
  static const unsigned LC_BOTH_EDGES = 2;

  const char *ReflectionBoundaryCondition::id = "ReflectionBc3D";

  ReflectionBoundaryCondition::ReflectionBoundaryCondition()
  {
  }

  void
  ReflectionBoundaryCondition::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis");
    else
      throw Lucee::Except("ReflectionBoundaryCondition::readInput: Must specify element to use using 'basis'");

    applyLeftEdge = applyRightEdge = false;
    std::string edgeStr = tbl.getString("edge");
    if (edgeStr == "lower")
      applyLeftEdge = true;
    else if (edgeStr == "upper")
      applyRightEdge = true;
    else if (edgeStr == "both")
      applyLeftEdge = applyRightEdge = true;
  }

  void
  ReflectionBoundaryCondition::initialize()
  {
    UpdaterIfc::initialize();

    unsigned nlocal = nodalBasis->getNumNodes();
    std::vector<unsigned> yRef(nlocal), xRef(nlocal);

    // reflection mapping
    rotMap.resize(nlocal);
    nodalBasis->getUpperReflectingBcMapping(0, yRef);
    nodalBasis->getLowerReflectingBcMapping(1, xRef);
    for (int i = 0; i < nlocal; i++)
      rotMap[i] = xRef[yRef[i]];
  }

  Lucee::UpdaterStatus
  ReflectionBoundaryCondition::update(double t)
  {
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    // Output distribution function
    Lucee::Field<3, double>& distf = this->getOut<Lucee::Field<3, double> >(0);

#ifdef HAVE_MPI
// This updater presently will not work in parallel. Eventually, we
// need to fix this, but it requires some MPI-fu to collect the values
// and put them on the appropriate processors. (Ammar Hakim,
// 8/06/2013)
    throw Lucee::Except("ReflectionBoundaryCondition does not work in parallel!");
#endif

    Lucee::Region<3, int> globalRgn = grid.getGlobalRegion();
    Lucee::FieldPtr<double> sknPtr = distf.createPtr(); // for skin-cell
    Lucee::FieldPtr<double> gstPtr = distf.createPtr(); // for ghost-cell

    unsigned nlocal = nodalBasis->getNumNodes();

    if (applyRightEdge == true)
    {
      int ix = globalRgn.getUpper(0)-1; // right skin cell x index

      // start with cell at maximum v_para, looping down to zero until desired flux is reached
      for (int js=globalRgn.getUpper(1)-1, jg=0; js>=0; --js, ++jg)
      {
        // Copy data into ghost after rotating skin cell data by 180 degrees
        for (int muIndex = globalRgn.getLower(2); muIndex < globalRgn.getUpper(2); muIndex++)
        {
          distf.setPtr(gstPtr, ix+1, jg, muIndex);
          distf.setPtr(sknPtr, ix, js, muIndex);

          for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
            gstPtr[componentIndex] = sknPtr[rotMap[componentIndex]];
        }
      }
    }

    if (applyLeftEdge == true)
    {
      int ix = globalRgn.getLower(0); // left skin cell x index

      // start with cell at maximum v_para, looping down to zero until desired flux is reached
      for (int js=globalRgn.getUpper(1)-1, jg=0; js>=0; --js, ++jg)
      {
        // Copy data into ghost after rotating skin cell data by 180 degrees
        for (int muIndex = globalRgn.getLower(2); muIndex < globalRgn.getUpper(2); muIndex++)
        {
          distf.setPtr(sknPtr, ix, js, muIndex);
          distf.setPtr(gstPtr, ix-1, jg, muIndex);

          for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
            gstPtr[componentIndex] = sknPtr[rotMap[componentIndex]];
        }
      }
    }

    return Lucee::UpdaterStatus();
  }

  void
  ReflectionBoundaryCondition::declareTypes()
  {
    // Distribution function output
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
  }
}
