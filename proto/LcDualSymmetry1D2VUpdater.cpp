/**
 * @file	LcDualSymmetry1D2VUpdater.cpp
 *
 * @brief	Simple updater that takes in a function W(y,py,px) and returns
 * and updater with the result W' = 0.5*[W(y,py,px) + W(y,-py,px)]
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDualSymmetry1D2VUpdater.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  const char *DualSymmetry1D2VUpdater::id = "DualSymmetry1D2VUpdater";

  DualSymmetry1D2VUpdater::DualSymmetry1D2VUpdater()
  {
  }

  DualSymmetry1D2VUpdater::~DualSymmetry1D2VUpdater()
  {
  }

  void
  DualSymmetry1D2VUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis");
    else
      throw Lucee::Except("DualSymmetry1D2VUpdater::readInput: Must specify element to use using 'basis'");
  }

  void
  DualSymmetry1D2VUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();
    
// get hold of grid
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();
    // local region to update
    Lucee::Region<3, int> localRgn = grid.getLocalRegion();

    unsigned nlocal = nodalBasis->getNumNodes();

    Lucee::RowMajorSequencer<3> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[3];
    seq.fillWithIndex(idx);
    nodalBasis->setIndex(idx);

    std::vector<unsigned> pyRef(nlocal), pxRef(nlocal);
    // Get reflection mapping after element has been reflected in py and px
    rotMap.resize(nlocal);
    nodalBasis->getUpperReflectingBcMapping(1, pyRef);
    nodalBasis->getUpperReflectingBcMapping(2, pxRef);
    for (int i = 0; i < nlocal; i++)
      rotMap[i] = pxRef[pyRef[i]];
  }

  Lucee::UpdaterStatus
  DualSymmetry1D2VUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    const Lucee::Field<3, double>& inpFld = this->getInp<Lucee::Field<3, double> >(0);
    Lucee::Field<3, double>& outFld = this->getOut<Lucee::Field<3, double> >(0);

    double dt = t-this->getCurrTime();
    
    Lucee::ConstFieldPtr<double> inpFldPtr  = inpFld.createConstPtr();
    Lucee::ConstFieldPtr<double> inpFldPtr_pair  = inpFld.createConstPtr();
    Lucee::FieldPtr<double> outFldPtr = outFld.createPtr();

    // Clear output field
    outFld = 0.0;

    // local region to index
    Lucee::Region<3, int> localRgn = grid.getLocalRegion();
    Lucee::Region<3, int> globalRgn = grid.getGlobalRegion();

    Lucee::RowMajorSequencer<3> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[3];
    double xc[3];
    seq.fillWithIndex(idx);
    nodalBasis->setIndex(idx);
    unsigned nlocal = nodalBasis->getNumNodes(); 

    // Loop over all position space cells in the local region
    for (int iy = localRgn.getLower(0); iy < localRgn.getUpper(0); iy++)
    {
      // Loop over all py and px. These must be on the same processor to avoid
      // dealing with MPI calls
      for (int ipy = globalRgn.getLower(1); ipy < globalRgn.getUpper(1); ipy++)
      {
        for (int ipx = globalRgn.getLower(2); ipx < globalRgn.getUpper(2); ipx++)
        {
          // Start with an arbitrary cell
          inpFld.setPtr(inpFldPtr, iy, ipy, ipx);
          // Find its pair according to W(y,-py,-px)
          inpFld.setPtr(inpFldPtr_pair, iy, globalRgn.getUpper(1) - 1 - ipy, 
              globalRgn.getUpper(2) - 1 - ipx);
          // Set output cell
          outFld.setPtr(outFldPtr, iy, ipy, ipx);

          for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
            outFldPtr[nodeIndex] = 0.5*(inpFldPtr[nodeIndex]+inpFldPtr_pair[ rotMap[nodeIndex] ]);
        }
      }
    }

    return Lucee::UpdaterStatus();
  }

  void
  DualSymmetry1D2VUpdater::declareTypes()
  {
    // Input field W(y,py,px)
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Output: W(y,py,px) = 0.5*( W(y,py,px) + W(y,-py,-px) )
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
  }
}
