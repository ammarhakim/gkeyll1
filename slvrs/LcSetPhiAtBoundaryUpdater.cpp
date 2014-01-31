/**
 * @file	LcSetPhiAtBoundaryUpdater.cpp
 *
 * @brief	Updater to modify an input phi(z) so that its value at the boundary is some value
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSetPhiAtBoundaryUpdater.h>
//#include <LcStructuredGridBase.h>
#include <LcMathPhysConstants.h>
// for cutoff velocities
#include <LcDynVector.h>

namespace Lucee
{
// set id for module system
  const char *SetPhiAtBoundaryUpdater::id = "SetPhiAtBoundaryUpdater";

  SetPhiAtBoundaryUpdater::SetPhiAtBoundaryUpdater()
    : UpdaterIfc()
  {
  }

  void 
  SetPhiAtBoundaryUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis");
    else
      throw Lucee::Except("SetPhiAtBoundaryUpdater::readInput: Must specify element to use using 'basis'");

    if (tbl.hasNumber("elcMass"))
      elcMass = tbl.getNumber("elcMass");
    else
      throw Lucee::Except("SetPhiAtBoundaryUpdater::readInput: Must specify elcMass");

    if (tbl.hasNumber("elcCharge"))
      elcCharge = tbl.getNumber("elcCharge");
    else
      throw Lucee::Except("SetPhiAtBoundaryUpdater::readInput: Must specify elcCharge");
  }

  void 
  SetPhiAtBoundaryUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus 
  SetPhiAtBoundaryUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<1>& grid
      = this->getGrid<Lucee::StructuredGridBase<1> >();

    const Lucee::Field<1, double>& phiIn = this->getInp<Lucee::Field<1, double> >(0);
    // Dynvector containing the cutoff velocities computed at one or both edges
    const Lucee::DynVector<double>& cutoffVIn = this->getInp<Lucee::DynVector<double> >(1);
    Lucee::Field<1, double>& phiOut = this->getOut<Lucee::Field<1, double> >(0);

    // Get exclusive node indices
    std::vector<int> exclusiveIndices;
    nodalBasis->getExclusiveNodeIndices(exclusiveIndices);
    //int nlocal = nodalBasis->getNumNodes();
    int nlocal = exclusiveIndices.size();

    Lucee::Region<1, int> extRegion = phiIn.getExtRegion();

    Lucee::ConstFieldPtr<double> phiInPtr = phiIn.createConstPtr();
    Lucee::FieldPtr<double> phiOutPtr = phiOut.createPtr();

    // Sheath potential (value of phi at right edge desired)
    double phiS;
    // Value of dPhi(x) on right boundary
    double dPhiRight;
    // Zero out output
    phiOutPtr = 0.0;

    // Get cutoff velocity on the right edge and calculate corresponding sheath potential
    std::vector<double> cutoffVelocities = cutoffVIn.getLastInsertedData();
    phiS = 0.5*elcMass*cutoffVelocities[1]*cutoffVelocities[1]/fabs(elcCharge);

    // Find offset that must be added to every cell
    // Assuming CG representation, this will be the ghost cell, zeroth component
    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion();
    phiIn.setPtr(phiInPtr, globalRgn.getUpper(0));
    dPhiRight = phiInPtr[0];

    // Loop over all cells
    for (int ix = extRegion.getLower(0); ix < extRegion.getUpper(0); ix++)
    {
      // Set inputs
      phiIn.setPtr(phiInPtr, ix);
      // Set outputs
      phiOut.setPtr(phiOutPtr, ix);

      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
        phiOutPtr[componentIndex] = phiInPtr[componentIndex] + phiS - dPhiRight;
    }
    
    return Lucee::UpdaterStatus();
  }

  void
  SetPhiAtBoundaryUpdater::declareTypes()
  {
    // takes two inputs: phi (CG field) and cutoff velocity vector
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::DynVector<double>));
    // returns one output: phi(x) with desired boundary conditions (CG field)
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }
}
