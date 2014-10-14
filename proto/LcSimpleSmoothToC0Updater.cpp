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

    Lucee::ConstFieldPtr<double> fInPtr = fIn.createConstPtr();
    Lucee::ConstFieldPtr<double> fInPtr_l = fIn.createConstPtr();
    Lucee::ConstFieldPtr<double> fInPtr_r = fIn.createConstPtr();
    Lucee::FieldPtr<double> fOutPtr = fOut.createPtr();
    Lucee::FieldPtr<double> fOutPtr_l = fOut.createPtr();
    Lucee::FieldPtr<double> fOutPtr_r = fOut.createPtr();

    fOut = 0.0;

    int idx[NDIM];
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);

    // Loop over edges in each direction
    for (int dir = 0; dir < NDIM; dir++)
    {
      // create sequencer to loop over *each* NDIM-1 slice in 'dir' direction
      Lucee::RowMajorSequencer<NDIM> seqLowerDim(localRgn.deflate(dir));
      // lower and upper bounds of NDIM-1 slice. (We need to make sure that flux
      // is computed for one edge outside domain interior)
      int sliceLower = localRgn.getLower(dir);
      int sliceUpper = localRgn.getUpper(dir)+1;

      int idxr[NDIM];
      int idxl[NDIM];

      // loop over each 1D slice
      while (seqLowerDim.step())
      {
        seqLowerDim.fillWithIndex(idxr);
        seqLowerDim.fillWithIndex(idxl);
        // loop over each edge
        for (int sliceIndex = sliceLower; sliceIndex < sliceUpper; sliceIndex++)
        { 
          idxr[dir] = sliceIndex;
          idxl[dir] = sliceIndex-1;

          fIn.setPtr(fInPtr_r, idxr);
          fIn.setPtr(fInPtr_l, idxl);

          fOut.setPtr(fOutPtr_r, idxr);
          fOut.setPtr(fOutPtr_l, idxl);
          // Loop over all nodes on a surface
          for (int nodeIndex = 0; nodeIndex < lowerNodeNums[dir].nums.size(); nodeIndex++)
          {
            int lowerNum = lowerNodeNums[dir].nums[nodeIndex];
            int upperNum = upperNodeNums[dir].nums[nodeIndex];

            fOutPtr_l[upperNum] = 0.25*(fInPtr_l[upperNum] + fInPtr_r[lowerNum]);
            fOutPtr_r[lowerNum] = 0.25*(fInPtr_l[upperNum] + fInPtr_r[lowerNum]);
          }
        }
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
  template class SimpleSmoothToC0Updater<1>;
  template class SimpleSmoothToC0Updater<2>;
  template class SimpleSmoothToC0Updater<3>;
  template class SimpleSmoothToC0Updater<4>;
  template class SimpleSmoothToC0Updater<5>;
}
