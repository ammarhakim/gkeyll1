/**
 * @file	LcASquaredProjectionUpdater.cpp
 *
 * @brief	Projects A(x)^2 onto the same basis functions that A(x) uses
 * Inputs and outputs assume DG fields, not CG fields.
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
  const char *ASquaredProjectionUpdater::id = "ASquaredProjectionUpdater";

  ASquaredProjectionUpdater::ASquaredProjectionUpdater()
    : UpdaterIfc()
  {
  }

  void 
  ASquaredProjectionUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis");
    else
      throw Lucee::Except("ASquaredProjectionUpdater::readInput: Must specify element to use using 'basis'");
  }

  void 
  ASquaredProjectionUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();
    // global region to update
    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion();

    Lucee::RowMajorSequencer<1> seq(globalRgn);
    seq.step(); // just to get to first index
    int idx[1];
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

    std::vector<double> gaussWeights(numQuadNodes);

    // Allocate Eigen matrices
    Eigen::MatrixXd interpMatrix(numQuadNodes, nlocal);
    Eigen::MatrixXd gaussOrdinates(numQuadNodes, 3);

    // Get the interpolation matrix for the volume quadrature points
    nodalBasis->getGaussQuadData(interpMatrixLucee, gaussOrdinatesLucee, gaussWeights);

    copyLuceeToEigen(massMatrixLucee, massMatrix);
    copyLuceeToEigen(interpMatrixLucee, interpMatrix);
    copyLuceeToEigen(gaussOrdinatesLucee, gaussOrdinates);

    // Compute and store the inverse of the mass matrix
    massMatrixInv = massMatrix.inverse();

    tripleProducts = std::vector<Eigen::MatrixXd>(nlocal);
    // Create and store the triple-product basis integrals

    for (int basisIndex = 0; basisIndex < nlocal; basisIndex++)
    {
      // Initialize matrices
      tripleProducts[basisIndex]  = Eigen::MatrixXd::Zero(nlocal, nlocal);
      for (int rowIndex = 0; rowIndex < nlocal; rowIndex++)
      {
        for (int colIndex = 0; colIndex < nlocal; colIndex++)
        {
          double integralResult = 0.0;

          for (int gaussNodeIndex = 0; gaussNodeIndex < numQuadNodes; gaussNodeIndex++)
          {
            integralResult += gaussWeights[gaussNodeIndex]*interpMatrix(gaussNodeIndex, basisIndex)*
              interpMatrix(gaussNodeIndex, rowIndex)*interpMatrix(gaussNodeIndex, colIndex);
          }

          tripleProducts[basisIndex](rowIndex, colIndex) = integralResult;
        }
      }
    }
  }

  Lucee::UpdaterStatus 
  ASquaredProjectionUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<1>& grid
      = this->getGrid<Lucee::StructuredGridBase<1> >();

    // Electron and ion densities
    const Lucee::Field<1, double>& aIn = this->getInp<Lucee::Field<1, double> >(0);
    Lucee::Field<1, double>& aSquaredOut = this->getOut<Lucee::Field<1, double> >(0);

    int nlocal = nodalBasis->getNumNodes();

    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion();

    Lucee::ConstFieldPtr<double> aPtr = aIn.createConstPtr();
    Lucee::FieldPtr<double> aSquaredPtr = aSquaredOut.createPtr();

    aSquaredPtr = 0.0;
    
    // Loop over all cells
    for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
    {
      // Set inputs
      aIn.setPtr(aPtr, ix);
      // Set outputs
      aSquaredOut.setPtr(aSquaredPtr, ix);

      Eigen::VectorXd aVec(nlocal);
      
      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
        aVec(componentIndex) = aPtr[componentIndex];

      Eigen::MatrixXd aProjections(nlocal, nlocal);

      for (int basisIndex = 0; basisIndex < nlocal; basisIndex++)
      {
        Eigen::VectorXd aProjectionsSingle = tripleProducts[basisIndex]*aVec;
        // Store components into rows of aProjections matrix
        for (int colIndex = 0; colIndex < nlocal; colIndex++)
          aProjections(basisIndex, colIndex) = aProjectionsSingle(colIndex);
      }

      // Multiply aProjections by aVec again to get a vector made up of
      // Int(phi_n*A(x)*A(x)) where n = component index, then multiply by
      // inverse of mass matrix to find the weights of A(x)^2
      Eigen::VectorXd aSquaredWeights = massMatrixInv*aProjections*aVec;

      // Copy into output pointer
      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
        aSquaredPtr[componentIndex] = aSquaredWeights(componentIndex);
    }

    return Lucee::UpdaterStatus();
  }

  void
  ASquaredProjectionUpdater::declareTypes()
  {
    // takes input A(x)
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // returns projection of A(x)^2 onto same order basis functions
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }

  void
  ASquaredProjectionUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}
