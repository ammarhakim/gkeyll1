/**
 * @file	LcMHDHamiltonianUpdater.cpp
 *
 * @brief	Computes the second-order gyrokinetic MHD Hamiltonian, without the mass factor
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcMHDHamiltonianUpdater.h>

namespace Lucee
{
// set id for module system
  const char *MHDHamiltonianUpdater::id = "MHDHamiltonianUpdater";

  MHDHamiltonianUpdater::MHDHamiltonianUpdater()
    : UpdaterIfc()
  {
  }

  void 
  MHDHamiltonianUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis");
    else
      throw Lucee::Except("MHDHamiltonianUpdater::readInput: Must specify element to use using 'basis'");
  }

  void 
  MHDHamiltonianUpdater::initialize()
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

    quadMassMatrix.resize(nlocal, nlocal, nlocal, nlocal);
    
    // Loop over nlocal x nlocal x nlocal x nlocal region
    int shape[4];
    int idx4d[4];
    blitz::TinyVector<int, 4> blitzCoord;

    for (int dimIndex = 0; dimIndex < 4; dimIndex++)
      shape[dimIndex] = nlocal;

    Lucee::Region<4, int> polyRegion(shape);
    Lucee::RowMajorSequencer<4> polySeq = RowMajorSequencer<4>(polyRegion);
    
    while(polySeq.step())
    {
      polySeq.fillWithIndex(idx4d);
      for (int dimIndex = 0; dimIndex < 4; dimIndex++)
        blitzCoord(dimIndex) = idx4d[dimIndex];

      double integralResult = 0.0;
      // Loop over each quadrature point, accumulating its contribution to integral
      for(int gaussNodeIndex = 0; gaussNodeIndex < numQuadNodes; gaussNodeIndex++)
      {
        double intermediateResult = gaussWeights[gaussNodeIndex];

        for(int termIndex = 0; termIndex < 4; termIndex++)
          intermediateResult *= interpMatrix(gaussNodeIndex, idx4d[termIndex]);

        integralResult += intermediateResult;
      }
      quadMassMatrix(blitzCoord) = integralResult;
    }

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
  MHDHamiltonianUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<1>& grid
      = this->getGrid<Lucee::StructuredGridBase<1> >();

    const Lucee::Field<1, double>& phiIn = this->getInp<Lucee::Field<1, double> >(0);
    const Lucee::Field<1, double>& nIn = this->getInp<Lucee::Field<1, double> >(1);
    Lucee::Field<1, double>& mhdHamiltonianOut = this->getOut<Lucee::Field<1, double> >(0);

    int nlocal = nodalBasis->getNumNodes();

    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion();

    Lucee::ConstFieldPtr<double> phiPtr = phiIn.createConstPtr();
    Lucee::ConstFieldPtr<double> nPtr = nIn.createConstPtr();
    Lucee::FieldPtr<double> mhdHamiltonianPtr = mhdHamiltonianOut.createPtr();

    mhdHamiltonianPtr = 0.0;
    
    // Loop over all cells
    for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
    {
      // Set inputs
      phiIn.setPtr(phiPtr, ix);
      nIn.setPtr(nPtr, ix);
      // Set outputs
      mhdHamiltonianOut.setPtr(mhdHamiltonianPtr, ix);

      Eigen::VectorXd phiVec(nlocal);
      Eigen::VectorXd nVec(nlocal);
      
      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
      {
        phiVec(componentIndex) = phiPtr[componentIndex];
        nVec(componentIndex) = nPtr[componentIndex];
      }

      // Construct LHS matrix
      Eigen::MatrixXd lhsProjections(nlocal, nlocal);
      for (int basisIndex = 0; basisIndex < nlocal; basisIndex++)
      {
        Eigen::VectorXd lhsProjectionsSingle = tripleProducts[basisIndex]*nVec;
        // Store components into rows of aProjections matrix
        for (int colIndex = 0; colIndex < nlocal; colIndex++)
          lhsProjections(basisIndex, colIndex) = lhsProjectionsSingle(colIndex);
      }

      // Construct RHS vector
      Eigen::VectorXd rhsVector(nlocal);
      for (int basisIndex = 0; basisIndex < nlocal; basisIndex++)
      {
        double integralResult = 0.0;
        for (int i = 0; i < nlocal; i++)
          for (int j = 0; j < nlocal; j++)
            for (int k = 0; k < nlocal; k++)
              integralResult += nVec(i)*phiVec(j)*phiVec(k)*quadMassMatrix(basisIndex, i, j, k);
        rhsVector(basisIndex) = integralResult;
      }

      Eigen::VectorXd mhdHamiltonianWeights = lhsProjections.inverse()*rhsVector;

      // Copy into output pointer
      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
        mhdHamiltonianPtr[componentIndex] = -0.5*mhdHamiltonianWeights(componentIndex);
    }

    return Lucee::UpdaterStatus();
  }

  void
  MHDHamiltonianUpdater::declareTypes()
  {
    // takes inputs phi1dDg, numDensity(x)
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // returns MHD hamiltonian
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }

  void
  MHDHamiltonianUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}
