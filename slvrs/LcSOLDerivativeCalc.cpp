/**
 * @file	LcSOLDerivativeCalc.cpp
 *
 * @brief	Computes the gradient of an input field in a particular direction
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLDerivativeCalc.h>

namespace Lucee
{
  template <> const char *SOLDerivativeCalc<1>::id = "SOLDerivativeCalc1D";
  template <> const char *SOLDerivativeCalc<2>::id = "SOLDerivativeCalc2D";
  template <> const char *SOLDerivativeCalc<3>::id = "SOLDerivativeCalc3D";
  template <> const char *SOLDerivativeCalc<4>::id = "SOLDerivativeCalc4D";
  template <> const char *SOLDerivativeCalc<5>::id = "SOLDerivativeCalc5D";

  template <unsigned NDIM>
  SOLDerivativeCalc<NDIM>::SOLDerivativeCalc()
  {
  }

  template <unsigned NDIM>
  void
  SOLDerivativeCalc<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("SOLDerivativeCalc::readInput: Must specify element to use using 'basis'");

    if (tbl.hasNumber("scaleFactor"))
      scaleFactor = tbl.getNumber("scaleFactor");
    else scaleFactor = 1.0;

    if (tbl.hasNumber("dir"))
      dir = (unsigned) tbl.getNumber("dir");
    else
      throw Lucee::Except("SOLDerivativeCalc::readInput: Must specify gradient direction using 'dir'");
  }

  template <unsigned NDIM>
  void
  SOLDerivativeCalc<NDIM>::initialize()
  {
    UpdaterIfc::initialize();

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
    grid.setIndex(idx);

    unsigned nlocal = nodalBasis->getNumNodes();
    // Store mass matrix inverse
    Lucee::Matrix<double> tempMatrix(nlocal, nlocal);
    nodalBasis->getMassMatrix(tempMatrix);
    Eigen::MatrixXd massMatrix(nlocal, nlocal);
    copyLuceeToEigen(tempMatrix, massMatrix);
    // Store stiffness matrix
    nodalBasis->getGradStiffnessMatrix(dir, tempMatrix);
    Eigen::MatrixXd gradStiffnessMatrix(nlocal, nlocal);
    copyLuceeToEigen(tempMatrix, gradStiffnessMatrix);
    // Compute and store differention matrix (net scale factor of 1/(0.5*grid.getDx(dir)))
    gradMatrix = massMatrix.inverse()*gradStiffnessMatrix.transpose();
    // Undo grid scale factors
    gradMatrix *= grid.getDx(dir);
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  SOLDerivativeCalc<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    // Input field
    const Lucee::Field<NDIM, double>& fieldIn = this->getInp<Lucee::Field<NDIM, double> >(0);
    // Output derivative field
    Lucee::Field<NDIM, double>& fieldOut = this->getOut<Lucee::Field<NDIM, double> >(0);
    
    fieldOut = 0.0; // clear out current contents

    Lucee::ConstFieldPtr<double> fieldInPtr = fieldIn.createConstPtr();
    Lucee::FieldPtr<double> fieldOutPtr = fieldOut.createPtr(); // Output pointer

    int idx[NDIM];
    unsigned nlocal = nodalBasis->getNumNodes();
    // Get ext region to loop over
    Lucee::Region<NDIM, int> localExtRgn = fieldOut.getExtRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localExtRgn);

    while (seq.step())
    {
      seq.fillWithIndex(idx);
      grid.setIndex(idx);
      fieldIn.setPtr(fieldInPtr, idx);
      fieldOut.setPtr(fieldOutPtr, idx);

      Eigen::VectorXd fieldInVec(nlocal);
      for (int i = 0; i < nlocal; i++)
        fieldInVec(i) = fieldInPtr[i];

      // Compute scale factor length(dir)
      double cellWidth = grid.getVolume()/grid.getSurfArea(dir);
      // Compute derivative
      Eigen::VectorXd fieldOutVec = gradMatrix*fieldInVec/cellWidth;

      // Write result
      for (int i = 0; i < nlocal; i++)
        fieldOutPtr[i] = scaleFactor*fieldOutVec(i);
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  SOLDerivativeCalc<NDIM>::declareTypes()
  {
    // input field
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    // output field with derivative applied in dir 'dir'
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void
  SOLDerivativeCalc<NDIM>::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
    {
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
      {
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
      }
    }
  }

// instantiations
  template class SOLDerivativeCalc<1>;
  template class SOLDerivativeCalc<2>;
  template class SOLDerivativeCalc<3>;
  template class SOLDerivativeCalc<4>;
  template class SOLDerivativeCalc<5>;
}
