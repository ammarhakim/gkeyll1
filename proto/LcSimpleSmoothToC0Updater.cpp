/**
 * @file	LcSimpleSmoothToC0Updater.cpp
 *
 * @brief	Updater to smooth a field by averaging nodal values on shared edges
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSimpleSmoothToC0Updater.h>

namespace Lucee
{
// set id for module system
  template <> const char *SimpleSmoothToC0Updater<1>::id = "SimpleSmoothToC01D";
  template <> const char *SimpleSmoothToC0Updater<2>::id = "SimpleSmoothToC02D";
  template <> const char *SimpleSmoothToC0Updater<3>::id = "SimpleSmoothToC03D";
  template <> const char *SimpleSmoothToC0Updater<4>::id = "SimpleSmoothToC04D";
  template <> const char *SimpleSmoothToC0Updater<5>::id = "SimpleSmoothToC05D";

  template <unsigned NDIM>
  SimpleSmoothToC0Updater<NDIM>::SimpleSmoothToC0Updater()
    : UpdaterIfc()
  {
  }
  
  template <unsigned NDIM>
  void 
  SimpleSmoothToC0Updater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("SimpleSmoothToC0Updater::readInput: Must specify element to use using 'basis'");
  }

  template <unsigned NDIM>
  void 
  SimpleSmoothToC0Updater<NDIM>::initialize()
  {
    // call base class method
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    // local region to update
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    
    // set index to first location in grid (this is okay as in this
    // updater we are assuming grid is uniform)
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    seq.step();
    int idx[NDIM];
    seq.fillWithIndex(idx);
    nodalBasis->setIndex(idx);

    // Get upper and lower node numbers
    for (int dir = 0; dir < NDIM; dir++)
    {
      lowerNodeNums[dir].nums.resize(nodalBasis->getNumSurfLowerNodes(dir));
      nodalBasis->getSurfLowerNodeNums(dir, lowerNodeNums[dir].nums);

      upperNodeNums[dir].nums.resize(nodalBasis->getNumSurfUpperNodes(dir));
      nodalBasis->getSurfUpperNodeNums(dir, upperNodeNums[dir].nums);
    }
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus 
  SimpleSmoothToC0Updater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    const Lucee::Field<NDIM, double>& fIn = this->getInp<Lucee::Field<NDIM, double> >(0);
    Lucee::Field<NDIM, double>& fOut = this->getOut<Lucee::Field<NDIM, double> >(0);

    // local region to update
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    Lucee::Region<NDIM, int> globalRgn = grid.getGlobalRegion();

    Lucee::ConstFieldPtr<double> fInPtr = fIn.createConstPtr();
    Lucee::ConstFieldPtr<double> fInPtr_l = fIn.createConstPtr();
    Lucee::ConstFieldPtr<double> fInPtr_r = fIn.createConstPtr();
    Lucee::FieldPtr<double> fOutPtr = fOut.createPtr();
    Lucee::FieldPtr<double> fOutPtr_l = fOut.createPtr();
    Lucee::FieldPtr<double> fOutPtr_r = fOut.createPtr();

    fOut = 0.0;

    int idx[NDIM];
    int idxl[NDIM];
    int idxr[NDIM];
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);

    // CODE ONLY WORKS WITH POLYORDER = 1 IN DIM = 2
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      seq.fillWithIndex(idxr);
      seq.fillWithIndex(idxl);

      fIn.setPtr(fInPtr, idx);
      fOut.setPtr(fOutPtr, idx);

      nodalBasis->setIndex(idx);

      // First initialize all nodes to 0.25 of current cell nodes
      for (int nodeIndex = 0; nodeIndex < nodalBasis->getNumNodes(); nodeIndex++)
        fOutPtr[nodeIndex] = 0.25*fInPtr[nodeIndex];

      // Loop over up,down,left,right neighboring cells
      for (int dir = 0; dir < NDIM; dir++)
      {
        // Get lower cell relative to this one in direction 'dir'
        idxl[dir] = idx[dir] - 1;
        // Get upper cell relative to this one in direction 'dir'
        idxr[dir] = idx[dir] + 1;

        fIn.setPtr(fInPtr_l, idxl);
        fIn.setPtr(fInPtr_r, idxr);

        // Loop over all nodes on a surface
        for (int nodeIndex = 0; nodeIndex < lowerNodeNums[dir].nums.size(); nodeIndex++)
        {
          int lowerNum = lowerNodeNums[dir].nums[nodeIndex];
          int upperNum = upperNodeNums[dir].nums[nodeIndex];

          fOutPtr[lowerNum] = fOutPtr[lowerNum] + 0.25*fInPtr_l[upperNum];
          fOutPtr[upperNum] = fOutPtr[upperNum] + 0.25*fInPtr_r[lowerNum];
        }

        // Restore indices to their original values before next iteration
        idxl[dir] = idx[dir];
        idxr[dir] = idx[dir];
      }

      // Loop over diagonal neighbors
      // Bottom left
      fIn.setPtr(fInPtr, idx[0]-1, idx[1]-1);
      fOutPtr[0] = fOutPtr[0] + 0.25*fInPtr[2];
      // Bottom right
      fIn.setPtr(fInPtr, idx[0]+1, idx[1]-1);
      fOutPtr[1] = fOutPtr[1] + 0.25*fInPtr[3];
      // Upper right
      fIn.setPtr(fInPtr, idx[0]+1, idx[1]+1);
      fOutPtr[2] = fOutPtr[2] + 0.25*fInPtr[0];
      // Upper left
      fIn.setPtr(fInPtr, idx[0]-1, idx[1]+1);
      fOutPtr[3] = fOutPtr[3] + 0.25*fInPtr[1];
    }

    return Lucee::UpdaterStatus();
  }
  
  template <unsigned NDIM>
  void
  SimpleSmoothToC0Updater<NDIM>::declareTypes()
  {
    // takes one input, a field to be smoothed
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    // returns one output, a smoothed field
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  // instantiations
  template class SimpleSmoothToC0Updater<1>;
  template class SimpleSmoothToC0Updater<2>;
  template class SimpleSmoothToC0Updater<3>;
  template class SimpleSmoothToC0Updater<4>;
  template class SimpleSmoothToC0Updater<5>;
}
