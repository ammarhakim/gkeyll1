/**
 * @file	LcNodalCopyFaceToInteriorUpdater.cpp
 *
 * @brief	Updater to copy 1D field to 2D field.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcNodalCopyFaceToInteriorUpdater.h>

namespace Lucee
{
  const char *NodalCopyFaceToInteriorUpdater::id = "NodalCopyFaceToInteriorUpdater";

  NodalCopyFaceToInteriorUpdater::NodalCopyFaceToInteriorUpdater()
    : Lucee::UpdaterIfc()
  {
  }

  void
  NodalCopyFaceToInteriorUpdater::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasNumber("dir"))
      dir = tbl.getNumber("dir");
    else
      throw Lucee::Except("NodalCopyFaceToInteriorUpdater::readInput: Must specify dir");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis1d"))
      nodalBasis1d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis1d");
    else
      throw Lucee::Except("NodalCopyFaceToInteriorUpdater::readInput: Must specify element to use using 'basis1d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis2d"))
      nodalBasis2d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis2d");
    else
      throw Lucee::Except("NodalCopyFaceToInteriorUpdater::readInput: Must specify element to use using 'basis2d'");

    // are common nodes shared?
    shareCommonNodes = false;
    if (tbl.hasBool("shareCommonNodes"))
      shareCommonNodes = tbl.getBool("shareCommonNodes");
  }

  void
  NodalCopyFaceToInteriorUpdater::initialize()
  {
    // call base class method
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();
    // Local region to update
    Lucee::Region<2, int> localRgn = grid.getLocalRegion();

    Lucee::RowMajorSequencer<2> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[1];
    seq.fillWithIndex(idx);
    nodalBasis2d->setIndex(idx);

    int nlocal = nodalBasis2d->getNumNodes();

    int nFaceNodes = nodalBasis2d->getNumSurfLowerNodes(dir);

    // Get the necessary mapping matrix
    Lucee::Matrix<double> mappingMatrixLucee(nlocal, nFaceNodes);
    nodalBasis2d->getLowerFaceToInteriorMapping(dir, mappingMatrixLucee);

    // Allocated Eigen matrices
    mappingMatrix = Eigen::MatrixXd(nlocal, nFaceNodes);

    // Copy Lucee matrix to Eigen matrix
    for (int rowIndex = 0; rowIndex < mappingMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < mappingMatrix.cols(); colIndex++)
        mappingMatrix(rowIndex, colIndex) = mappingMatrixLucee(rowIndex, colIndex);

    // If shareCommonNodes == true, we need to do additional work
    if (shareCommonNodes == true)
    {
      // Get exclusive node indices
      std::vector<int> exclusiveIndices;
      nodalBasis2d->getExclusiveNodeIndices(exclusiveIndices);

      // Temporary matrix to hold mapping
      Eigen::MatrixXd tempMappingMatrix(exclusiveIndices.size(), nFaceNodes);

      // Based on the nodes listed in getExclusiveNodeIndices, only
      // keep these rows from the mappingMatrix
      for (int rowIndex = 0; rowIndex < exclusiveIndices.size(); rowIndex++)
      {
        // Copy all elements in row = exclusiveIndices[rowIndex] into row = rowIndex
        // of the temporary matrix
        for (int colIndex = 0; colIndex < mappingMatrix.cols(); colIndex++)
          tempMappingMatrix(rowIndex, colIndex) = mappingMatrix(exclusiveIndices[rowIndex], colIndex);
      }

      mappingMatrix = tempMappingMatrix;
    }
  }

  Lucee::UpdaterStatus
  NodalCopyFaceToInteriorUpdater::update(double t)
  {
    // get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();

    // get input field (1d)
    const Lucee::Field<1, double>& fld1d = this->getInp<Lucee::Field<1, double> >(0);
    // get output field (2D)
    Lucee::Field<2, double>& fld2d = this->getOut<Lucee::Field<2, double> >(0);

    // local region to update (This is the 2D region. The 1D region is
    // assumed to have the same cell layout as the X-direction of the 2D region)
    Lucee::Region<2, int> localRgn = grid.getLocalRegion();

    Lucee::ConstFieldPtr<double> fld1dPtr = fld1d.createConstPtr();
    Lucee::FieldPtr<double> fld2dPtr = fld2d.createPtr();

    int nFaceNodes = nodalBasis2d->getNumSurfLowerNodes(dir);

    // loop over all X-direction cells
    for (int i=localRgn.getLower(0)-1; i<localRgn.getUpper(0)+1; ++i)
    {
      fld1d.setPtr(fld1dPtr, i);
      Eigen::VectorXd vectorFace(mappingMatrix.cols());

      if (shareCommonNodes == true)
      {
        // Extract from field at this cell
        std::vector<double> completeFaceData(nFaceNodes);
        nodalBasis1d->setIndex(i);
        nodalBasis1d->extractFromField(fld1d, completeFaceData);
        
        // Copy into Eigen vectors to do matrix multiplies
        for (int componentIndex = 0; componentIndex < mappingMatrix.cols(); componentIndex++)
          vectorFace(componentIndex) = completeFaceData[componentIndex];
      }
      else
      {
        // No extract from field needed. Already have complete face data
        // Copy into Eigen vectors to do matrix multiplies
        for( int componentIndex = 0; componentIndex < mappingMatrix.cols(); componentIndex++)
          vectorFace(componentIndex) = fld1dPtr[componentIndex];
      }
      // Apply this mapping to get 2-D nodal data
      Eigen::VectorXd vectorVolume = mappingMatrix*vectorFace;

      for (int j=localRgn.getLower(1)-1; j<localRgn.getUpper(1)+1; ++j)
      {
        fld2d.setPtr(fld2dPtr, i, j);
        for (int componentIndex = 0; componentIndex < vectorVolume.rows(); componentIndex++)
          fld2dPtr[componentIndex] = vectorVolume(componentIndex);
      }
    }

    return Lucee::UpdaterStatus();
  }

  void
  NodalCopyFaceToInteriorUpdater::declareTypes()
  {
    // Input: Face field
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // Output: Volume field
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }
}
