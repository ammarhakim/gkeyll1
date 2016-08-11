/**
 * @file	LcASquaredProjectionUpdater.cpp
 *
 * @brief	Updater to project a product of two fields onto the original basis set using quadrature
 * Inputs and outputs assume DG fields, not CG fields.
 * 6-19-16: I rewrote this updater, which originally computed A^2 to compute A*B as a way to make it more
 * useful. This will break existing Lua scripts that make use of it, but the changes will be minor. (ELS)
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcASquaredProjectionUpdater.h>

namespace Lucee
{
// set id for module system
  template <> const char *ASquaredProjectionUpdater<1>::id = "ASquaredProjection1D";
  template <> const char *ASquaredProjectionUpdater<2>::id = "ASquaredProjection2D";
  template <> const char *ASquaredProjectionUpdater<3>::id = "ASquaredProjection3D";

  template <unsigned NDIM>
  ASquaredProjectionUpdater<NDIM>::ASquaredProjectionUpdater()
    : UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void 
  ASquaredProjectionUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("ASquaredProjectionUpdater::readInput: Must specify element to use using 'basis'");
  }

  template <unsigned NDIM>
  void 
  ASquaredProjectionUpdater<NDIM>::initialize()
  {
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    // local region to update
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[NDIM];
    seq.fillWithIndex(idx);
    nodalBasis->setIndex(idx);
    
    int nlocal = nodalBasis->getNumNodes();
    int numQuadNodes = nodalBasis->getNumGaussNodes();

    // Get mass matrix and then copy to Eigen format
    Lucee::Matrix<double> massMatrixLucee(nlocal, nlocal);

    Eigen::MatrixXd massMatrix(nlocal, nlocal);
    
    nodalBasis->getMassMatrix(massMatrixLucee);
    
    // Get interpolation matrix, gaussian quadrature points, and weights
    Lucee::Matrix<double> interpMatrixLucee(numQuadNodes, nlocal);

    Lucee::Matrix<double> gaussOrdinatesLucee(numQuadNodes, 3);

    gaussWeights = std::vector<double>(numQuadNodes);

    // Allocate Eigen matrices
    interpMatrix = Eigen::MatrixXd(numQuadNodes, nlocal);
    Eigen::MatrixXd gaussOrdinates(numQuadNodes, 3);

    // Get the interpolation matrix for the volume quadrature points
    nodalBasis->getGaussQuadData(interpMatrixLucee, gaussOrdinatesLucee, gaussWeights);

    copyLuceeToEigen(massMatrixLucee, massMatrix);
    copyLuceeToEigen(interpMatrixLucee, interpMatrix);
    copyLuceeToEigen(gaussOrdinatesLucee, gaussOrdinates);

    // Compute and store the inverse of the mass matrix
    massMatrixInv = massMatrix.inverse();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus 
  ASquaredProjectionUpdater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    // Electron and ion densities
    const Lucee::Field<NDIM, double>& aIn = this->getInp<Lucee::Field<NDIM, double> >(0);
    const Lucee::Field<NDIM, double>& bIn = this->getInp<Lucee::Field<NDIM, double> >(1);
    Lucee::Field<NDIM, double>& productOut = this->getOut<Lucee::Field<NDIM, double> >(0);

    int nlocal = nodalBasis->getNumNodes();

    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    Lucee::ConstFieldPtr<double> aPtr = aIn.createConstPtr();
    Lucee::ConstFieldPtr<double> bPtr = bIn.createConstPtr();
    Lucee::FieldPtr<double> productPtr = productOut.createPtr();

    productPtr = 0.0;
    
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    int idx[NDIM];
    Eigen::VectorXd aVec(nlocal);
    Eigen::VectorXd bVec(nlocal);
    Eigen::VectorXd aAtQuad(gaussWeights.size());
    Eigen::VectorXd bAtQuad(gaussWeights.size());
    Eigen::VectorXd projectionVec(nlocal);
    Eigen::VectorXd solutionVec(nlocal);

    // Loop over all local cells
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      // Set inputs
      aIn.setPtr(aPtr, idx);
      bIn.setPtr(bPtr, idx);
      // Set outputs
      productOut.setPtr(productPtr, idx);
      
      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
      {
        aVec(componentIndex) = aPtr[componentIndex];
        bVec(componentIndex) = bPtr[componentIndex];
      }

      aAtQuad = interpMatrix*aVec;
      bAtQuad = interpMatrix*bVec;

      // Compute projection of product A*B onto each basis function
      for (int basisIndex = 0; basisIndex < nlocal; basisIndex++)
      {
        double integralResult = 0.0;
        for (int quadIndex = 0; quadIndex < gaussWeights.size(); quadIndex++)
          integralResult += gaussWeights[quadIndex]*interpMatrix(quadIndex,basisIndex)*
            aAtQuad(quadIndex)*bAtQuad(quadIndex);
        projectionVec(basisIndex) = integralResult;
      }
      solutionVec = massMatrixInv*projectionVec;

      // Copy into output pointer
      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
        productPtr[componentIndex] = solutionVec(componentIndex);
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  ASquaredProjectionUpdater<NDIM>::declareTypes()
  {
    // takes input A(x), B(x)
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    // returns projection of A(x)*B(x) onto same order basis functions
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void
  ASquaredProjectionUpdater<NDIM>::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }

  // instantiations
  template class ASquaredProjectionUpdater<1>;
  template class ASquaredProjectionUpdater<2>;
  template class ASquaredProjectionUpdater<3>;
}
