/**
 * @file	LcPositivityUpdater.cpp
 *
 * @brief	Updater to enforce positivity preservation.
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
#include <LcPositivityUpdater.h>
#include <LcStructuredGridBase.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
// set id for module system
  const char *PositivityUpdater::id = "PositivityUpdater";

  PositivityUpdater::PositivityUpdater()
    : UpdaterIfc()
  {
  }

  void 
  PositivityUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis");
    else
      throw Lucee::Except("PositivityUpdater::readInput: Must specify element to use using 'basis'");
  }

  void 
  PositivityUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();

    //if (nodalBasis->getNumNodes() > 4)
    //  throw Lucee::Except("PositivityUpdater::initialize: Only polyOrder = 1 is supported right now.");
  }

  Lucee::UpdaterStatus 
  PositivityUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<2>& grid
      = this->getGrid<Lucee::StructuredGridBase<2> >();

    // Input field that can have negative parts
    const Lucee::Field<2, double>& distfIn = this->getInp<Lucee::Field<2, double> >(0);
    // Output field that has been modified
    Lucee::Field<2, double>& distfOut = this->getOut<Lucee::Field<2, double> >(0);

    int nlocal = nodalBasis->getNumNodes();

    Lucee::ConstFieldPtr<double> distfPtr = distfIn.createConstPtr();
    Lucee::FieldPtr<double> distfModifiedPtr = distfOut.createPtr();
    
    //Lucee::Region<2, int> globalRgn = grid.getGlobalRegion();
    Lucee::Region<2, int> extRegion = distfIn.getExtRegion();

    distfModifiedPtr = 0.0;

    for (int ix = extRegion.getLower(0); ix < extRegion.getUpper(0); ix++)
    {
      for (int iv = extRegion.getLower(1); iv < extRegion.getUpper(1); iv++)
      {
        // Set inputs
        distfIn.setPtr(distfPtr, ix, iv);
        // Set outputs
        distfOut.setPtr(distfModifiedPtr, ix, iv);

        // Check to make sure modifying values will maintain positivity...
        /*double cellAvg = 0.5*(distfPtr[0] + distfPtr[3]);
        if (cellAvg < 0.0)
          throw Lucee::Except("PositivityUpdater::update: cellAvg is negative...");

        cellAvg = 0.5*(distfPtr[1] + distfPtr[2]);
        if (cellAvg < 0.0)
          throw Lucee::Except("PositivityUpdater::update: cellAvg is negative...");*/
        /*for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
        {
          distfModifiedPtr[componentIndex] = distfPtr[componentIndex];
          if (distfPtr[componentIndex] < 0.0)
            distfModifiedPtr[componentIndex] = 0.0;
        }*/

        if (distfPtr[3] + distfPtr[0] < 0.0)
        {
          distfModifiedPtr[0] = 0.0;
          distfModifiedPtr[3] = 0.0;
        }
        else if (distfPtr[0] < 0.0 && distfPtr[3] > 0.0)
        {
          distfModifiedPtr[0] = 0.0;
          distfModifiedPtr[3] = distfPtr[3] + distfPtr[0];
        }
        else if (distfPtr[3] < 0.0 && distfPtr[0] > 0.0)
        {
          distfModifiedPtr[0] = distfPtr[3] + distfPtr[0];
          distfModifiedPtr[3] = 0.0;
        }
        else
        {
          distfModifiedPtr[0] = distfPtr[0];
          distfModifiedPtr[3] = distfPtr[3];
        }

        if (distfPtr[1] + distfPtr[2] < 0.0)
        {
          distfModifiedPtr[1] = 0.0;
          distfModifiedPtr[2] = 0.0;
        }
        else if (distfPtr[1] < 0.0 && distfPtr[2] > 0.0)
        {
          distfModifiedPtr[1] = 0.0;
          distfModifiedPtr[2] = distfPtr[1] + distfPtr[2];
        }
        else if (distfPtr[2] < 0.0 && distfPtr[1] > 0.0)
        {
          distfModifiedPtr[1] = distfPtr[1] + distfPtr[2];
          distfModifiedPtr[2] = 0.0;
        }
        else
        {
          distfModifiedPtr[1] = distfPtr[1];
          distfModifiedPtr[2] = distfPtr[2];
        }
      }
    }

    return Lucee::UpdaterStatus();
  }

  void
  PositivityUpdater::declareTypes()
  {
    // takes one input: original field
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    // returns one output: modified field (e.g. distribution function)
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }
}
