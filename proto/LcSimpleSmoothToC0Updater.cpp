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
  template <> const char *SimpleSmoothToC0Updater<2>::id = "SimpleSmoothToC02D";
  template <> const char *SimpleSmoothToC0Updater<3>::id = "SimpleSmoothToC03D";

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

    if (tbl.hasNumber("polyOrder"))
      polyOrder = tbl.getNumber("polyOrder");
    else
      throw Lucee::Except("SimpleSmoothToC0Updater::readInput: Must specify polyOrder");

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

      if (NDIM == 2)
      {
        // First initialize all nodes to 0.5 of current cell nodes
        for (int nodeIndex = 0; nodeIndex < nodalBasis->getNumNodes(); nodeIndex++)
          fOutPtr[nodeIndex] = 0.5*fInPtr[nodeIndex];

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

            fOutPtr[lowerNum] = fOutPtr[lowerNum] + 0.5*fInPtr_l[upperNum];
            fOutPtr[upperNum] = fOutPtr[upperNum] + 0.5*fInPtr_r[lowerNum];
          }

          // Restore indices to their original values before next iteration
          idxl[dir] = idx[dir];
          idxr[dir] = idx[dir];
        }

        // Loop over diagonal neighbors
        if (polyOrder == 1)
        {
          // Bottom left
          fIn.setPtr(fInPtr, idx[0]-1, idx[1]-1);
          fOutPtr[0] = 0.5*fOutPtr[0] + 0.25*fInPtr[3];
          // Bottom right
          fIn.setPtr(fInPtr, idx[0]+1, idx[1]-1);
          fOutPtr[1] = 0.5*fOutPtr[1] + 0.25*fInPtr[2];
          // Upper left
          fIn.setPtr(fInPtr, idx[0]-1, idx[1]+1);
          fOutPtr[2] = 0.5*fOutPtr[2] + 0.25*fInPtr[1]; 
          // Upper right
          fIn.setPtr(fInPtr, idx[0]+1, idx[1]+1);
          fOutPtr[3] = 0.5*fOutPtr[3] + 0.25*fInPtr[0];
        }
        else if (polyOrder == 2)
        {
          // Bottom left
          fIn.setPtr(fInPtr, idx[0]-1, idx[1]-1);
          fOutPtr[0] = 0.5*fOutPtr[0] + 0.25*fInPtr[7];
          // Bottom right
          fIn.setPtr(fInPtr, idx[0]+1, idx[1]-1);
          fOutPtr[2] = 0.5*fOutPtr[2] + 0.25*fInPtr[5];
          // Upper left
          fIn.setPtr(fInPtr, idx[0]-1, idx[1]+1);
          fOutPtr[5] = 0.5*fOutPtr[5] + 0.25*fInPtr[2];
          // Upper right
          fIn.setPtr(fInPtr, idx[0]+1, idx[1]+1);
          fOutPtr[7] = 0.5*fOutPtr[7] + 0.25*fInPtr[0];
        }
      }
      else if (NDIM == 3)
      {
        // Three types of loops needed: shared faces, shared edges, and shared vertices
        
        // Each corner node value is averaged with its eight neighbors
        double scaleFactor = 1.0/8.0;

        // First initialize all nodes to 1/8 of current cell nodes
        for (int nodeIndex = 0; nodeIndex < nodalBasis->getNumNodes(); nodeIndex++)
          fOutPtr[nodeIndex] = scaleFactor*fInPtr[nodeIndex];

        // Loop over six cells that share a common face
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

            fOutPtr[lowerNum] = fOutPtr[lowerNum] + scaleFactor*fInPtr_l[upperNum];
            fOutPtr[upperNum] = fOutPtr[upperNum] + scaleFactor*fInPtr_r[lowerNum];
          }

          // Restore indices to their original values before next iteration
          idxl[dir] = idx[dir];
          idxr[dir] = idx[dir];
        }

        // Loop over twelve cells that share a common edge
        fIn.setPtr(fInPtr, idx[0]-1, idx[1]-1, idx[2]);
        fOutPtr[0] = fOutPtr[0] + scaleFactor*fInPtr[3];
        fOutPtr[4] = fOutPtr[4] + scaleFactor*fInPtr[7];
        fIn.setPtr(fInPtr, idx[0]+1, idx[1]-1, idx[2]);
        fOutPtr[1] = fOutPtr[1] + scaleFactor*fInPtr[2];
        fOutPtr[5] = fOutPtr[5] + scaleFactor*fInPtr[6];
        fIn.setPtr(fInPtr, idx[0]-1, idx[1]+1, idx[2]);
        fOutPtr[2] = fOutPtr[2] + scaleFactor*fInPtr[1];
        fOutPtr[6] = fOutPtr[6] + scaleFactor*fInPtr[5];
        fIn.setPtr(fInPtr, idx[0]+1, idx[1]+1, idx[2]);
        fOutPtr[3] = fOutPtr[3] + scaleFactor*fInPtr[0];
        fOutPtr[7] = fOutPtr[7] + scaleFactor*fInPtr[4];

        fIn.setPtr(fInPtr, idx[0]-1, idx[1], idx[2]-1);
        fOutPtr[0] = fOutPtr[0] + scaleFactor*fInPtr[5];
        fOutPtr[2] = fOutPtr[2] + scaleFactor*fInPtr[7];
        fIn.setPtr(fInPtr, idx[0]+1, idx[1], idx[2]-1);
        fOutPtr[1] = fOutPtr[1] + scaleFactor*fInPtr[4];
        fOutPtr[3] = fOutPtr[3] + scaleFactor*fInPtr[6];
        fIn.setPtr(fInPtr, idx[0]-1, idx[1], idx[2]+1);
        fOutPtr[4] = fOutPtr[4] + scaleFactor*fInPtr[1];
        fOutPtr[6] = fOutPtr[6] + scaleFactor*fInPtr[3];
        fIn.setPtr(fInPtr, idx[0]+1, idx[1], idx[2]+1);
        fOutPtr[5] = fOutPtr[5] + scaleFactor*fInPtr[0];
        fOutPtr[7] = fOutPtr[7] + scaleFactor*fInPtr[2];

        fIn.setPtr(fInPtr, idx[0], idx[1]-1, idx[2]-1);
        fOutPtr[0] = fOutPtr[0] + scaleFactor*fInPtr[6];
        fOutPtr[1] = fOutPtr[1] + scaleFactor*fInPtr[7];
        fIn.setPtr(fInPtr, idx[0], idx[1]+1, idx[2]-1);
        fOutPtr[2] = fOutPtr[2] + scaleFactor*fInPtr[4];
        fOutPtr[3] = fOutPtr[3] + scaleFactor*fInPtr[5];
        fIn.setPtr(fInPtr, idx[0], idx[1]-1, idx[2]+1);
        fOutPtr[4] = fOutPtr[4] + scaleFactor*fInPtr[2];
        fOutPtr[5] = fOutPtr[5] + scaleFactor*fInPtr[3];
        fIn.setPtr(fInPtr, idx[0], idx[1]+1, idx[2]+1);
        fOutPtr[6] = fOutPtr[6] + scaleFactor*fInPtr[0];
        fOutPtr[7] = fOutPtr[7] + scaleFactor*fInPtr[1];

        // Loop over eight cells that only share a corner intersection
        fIn.setPtr(fInPtr, idx[0]-1, idx[1]-1, idx[2]-1);
        fOutPtr[0] = fOutPtr[0] + scaleFactor*fInPtr[7];
        fIn.setPtr(fInPtr, idx[0]+1, idx[1]-1, idx[2]-1);
        fOutPtr[1] = fOutPtr[1] + scaleFactor*fInPtr[6];
        fIn.setPtr(fInPtr, idx[0]-1, idx[1]+1, idx[2]-1);
        fOutPtr[2] = fOutPtr[2] + scaleFactor*fInPtr[5];
        fIn.setPtr(fInPtr, idx[0]+1, idx[1]+1, idx[2]-1);
        fOutPtr[3] = fOutPtr[3] + scaleFactor*fInPtr[4];
        fIn.setPtr(fInPtr, idx[0]-1, idx[1]-1, idx[2]+1);
        fOutPtr[4] = fOutPtr[4] + scaleFactor*fInPtr[3];
        fIn.setPtr(fInPtr, idx[0]+1, idx[1]-1, idx[2]+1);
        fOutPtr[5] = fOutPtr[5] + scaleFactor*fInPtr[2];
        fIn.setPtr(fInPtr, idx[0]-1, idx[1]+1, idx[2]+1);
        fOutPtr[6] = fOutPtr[6] + scaleFactor*fInPtr[1];
        fIn.setPtr(fInPtr, idx[0]+1, idx[1]+1, idx[2]+1);
        fOutPtr[7] = fOutPtr[7] + scaleFactor*fInPtr[0];
      }
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
  template class SimpleSmoothToC0Updater<2>;
  template class SimpleSmoothToC0Updater<3>;
}
