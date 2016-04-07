/**
 * @file	LcSOLSetPotentialAtBoundary.cpp
 *
 * @brief	Sets Dirchlet boundary conditions on potential for 5D SOL simulations
 * Given the 2d fields on the lower and upper z surfaces, puts them into a 3d field
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLSetPotentialAtBoundary.h>
#include <LcGlobals.h>
#include <LcLinAlgebra.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <cmath>
#include <vector>

namespace Lucee
{
  const char *SOLSetPotentialAtBoundary::id = "SOLSetPotentialAtBoundary";

  SOLSetPotentialAtBoundary::SOLSetPotentialAtBoundary()
    : Lucee::UpdaterIfc()
  {
  }

  SOLSetPotentialAtBoundary::~SOLSetPotentialAtBoundary()
  {
  }
  
  void
  SOLSetPotentialAtBoundary::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    // get hold of 3D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis");
    else
      throw Lucee::Except("SOLSetPotentialAtBoundary::readInput: Must specify 3D element to use using 'basis'");

    // by default, only set potential to a constant on the upper x edge
    // here is an option to also set potential to a constant lower x edge
    if (tbl.hasBool("applyLowerEdge"))
      applyLowerEdge = tbl.getBool("applyLowerEdge");
    else applyLowerEdge = false;
  }

  void
  SOLSetPotentialAtBoundary::initialize()
  {
    // call base class method
    Lucee::UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus
  SOLSetPotentialAtBoundary::update(double t)
  {
    
    // get hold of grid
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    // input 2d potential on lower z surface
    const Lucee::Field<2, double>& phiZLowerIn = this->getInp<Lucee::Field<2, double> >(0);
    // input 2d potential on upper z surface
    const Lucee::Field<2, double>& phiZUpperIn = this->getInp<Lucee::Field<2, double> >(1);
    // input DynVector containing value to set phi on upper 0 surface
    const Lucee::DynVector<double>& phiXUpperIn = this->getInp<Lucee::DynVector<double> >(2);
    // output 3d potential
    Lucee::Field<3, double>& phi3dOut = this->getOut<Lucee::Field<3, double> >(0);

    // Local region to update
    Lucee::Region<3, int> localRgn = grid.getLocalRegion();
    // Global region to tell if we are on domain boundary
    Lucee::Region<3, int> globalRgn = grid.getGlobalRegion();

    // iterators into fields
    Lucee::ConstFieldPtr<double> phiZLowerInPtr = phiZLowerIn.createConstPtr();
    Lucee::ConstFieldPtr<double> phiZUpperInPtr = phiZUpperIn.createConstPtr();
    Lucee::FieldPtr<double> phi3dOutPtr = phi3dOut.createPtr();
 
    int idx[3];
    Lucee::RowMajorSequencer<3> seq(localRgn);
    unsigned nlocal3d = nodalBasis3d->getNumNodes();

    // Get value from phiXUpperIn
    std::vector<double> phiXUpperVec = phiXUpperIn.getLastInsertedData();
    //double phiXUpperVal = phiXUpperVec[0];

    // Loop over every cell in the local region
    // Only set a value if it is on a Dirchlet BC surface
    while(seq.step())
    {
      seq.fillWithIndex(idx);
      grid.setIndex(idx);
      phi3dOut.setPtr(phi3dOutPtr, idx);

      // Set value on upper 0 surface to dynvector value
      if ( idx[0] == globalRgn.getUpper(0)-1 )
      {
        std::vector<int> nodeNums(nodalBasis3d->getNumSurfUpperNodes(0));
        nodalBasis3d->getSurfUpperNodeNums(0, nodeNums);

        for (int i = 0; i < nodeNums.size(); i++)
          phi3dOutPtr[nodeNums[i]] = phiXUpperVec[0];
      }

      if (applyLowerEdge == true && idx[0] == globalRgn.getLower(0))
      {
        std::vector<int> nodeNums(nodalBasis3d->getNumSurfLowerNodes(0));
        nodalBasis3d->getSurfLowerNodeNums(0, nodeNums);

        for (int i = 0; i < nodeNums.size(); i++)
          phi3dOutPtr[nodeNums[i]] = phiXUpperVec[1];
      }

      // Set value on lower 2 surface to phiZLower values
      if ( idx[2] == globalRgn.getLower(2) )
      {
        phiZLowerIn.setPtr(phiZLowerInPtr, idx[0], idx[1]);
        std::vector<int> nodeNums(nodalBasis3d->getNumSurfLowerNodes(2));
        nodalBasis3d->getSurfLowerNodeNums(2, nodeNums);

        for (int i = 0; i < nodeNums.size(); i++)
          phi3dOutPtr[nodeNums[i]] = phiZLowerInPtr[i];
      }

      // Set value on upper 2 surface to phiZUpper values
      if ( idx[2] == globalRgn.getUpper(2)-1 )
      {
        phiZUpperIn.setPtr(phiZUpperInPtr, idx[0], idx[1]);
        std::vector<int> nodeNums(nodalBasis3d->getNumSurfUpperNodes(2));
        nodalBasis3d->getSurfUpperNodeNums(2, nodeNums);

        for (int i = 0; i < nodeNums.size(); i++)
          phi3dOutPtr[nodeNums[i]] = phiZUpperInPtr[i];
      }
    }
    
    return Lucee::UpdaterStatus();
  }

  void
  SOLSetPotentialAtBoundary::declareTypes()
  {
    // Input 2d potential on lower z surface
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    // Input 2d potential on upper z surface
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    // Input DynVector containing value for phi on upper x surface
    this->appendInpVarType(typeid(Lucee::DynVector<double>));
    // Output 3d potential containing Dirchlet boundary values
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
  }
}
