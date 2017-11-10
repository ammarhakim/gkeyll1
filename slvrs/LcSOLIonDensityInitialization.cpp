/**
 * @file	LcSOLIonDensityInitialization.cpp
 *
 * @brief	Updater to compute an ion density given an electron density in accordance
 * with force balance
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLIonDensityInitialization.h>
#include <LcMathPhysConstants.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <vector>

namespace Lucee
{
// set id for module system
  const char *SOLIonDensityInitialization::id = "SOLIonDensityInitialization";
  const unsigned NDIM = 1;

  SOLIonDensityInitialization::SOLIonDensityInitialization()
    : UpdaterIfc()
  {
  }

  void 
  SOLIonDensityInitialization::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("SOLIonDensityInitialization::readInput: Must specify element to use using 'basis'");

    if (tbl.hasNumber("kPerpTimesRho"))
      kPerpTimesRho = tbl.getNumber("kPerpTimesRho");
    else
      throw Lucee::Except("SOLIonDensityInitialization::readInput: Must specify kPerpTimesRho");

    if (tbl.hasNumber("Te0"))
      Te0 = tbl.getNumber("Te0");
  }

  void 
  SOLIonDensityInitialization::initialize()
  {
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[NDIM];
    seq.fillWithIndex(idx);
    nodalBasis->setIndex(idx);
    
    int nlocal = nodalBasis->getNumNodes();
    int numQuadNodes = nodalBasis->getNumGaussNodes();

    // Get interpolation matrix, gaussian quadrature points, and weights
    Lucee::Matrix<double> interpMatrixLucee(numQuadNodes, nlocal);
    Lucee::Matrix<double> gaussOrdinatesLucee(numQuadNodes, 3);
    gaussWeights = std::vector<double>(numQuadNodes);
    // Get the interpolation matrix for the volume quadrature points
    nodalBasis->getGaussQuadData(interpMatrixLucee, gaussOrdinatesLucee, gaussWeights);

    // Allocate Eigen matrices
    interpMatrix = Eigen::MatrixXd(numQuadNodes, nlocal);

    copyLuceeToEigen(interpMatrixLucee, interpMatrix);
  }

  Lucee::UpdaterStatus 
  SOLIonDensityInitialization::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    // Electron and ion densities
    const Lucee::Field<NDIM, double>& nElcIn = this->getInp<Lucee::Field<NDIM, double> >(0);
    Lucee::Field<NDIM, double>& nIonOut = this->getOut<Lucee::Field<NDIM, double> >(0);

    int nlocal = nodalBasis->getNumNodes();

    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    Lucee::ConstFieldPtr<double> nElcInPtr = nElcIn.createConstPtr();
    Lucee::FieldPtr<double> nIonOutPtr = nIonOut.createPtr();

    nIonOutPtr = 0.0;

    int numCells = localRgn.getUpper(0)-localRgn.getLower(0);
    double domainLength = numCells*grid.getDx(0);

    Eigen::VectorXd nIonPrev = Eigen::VectorXd::Zero(numCells*nlocal);

    // Calculate integrated electron density and electron-density-weighted ln(n_e) (one time)
    double intElcDensityLocal = 0.0;
    double intElcDensity = 0.0;
    double avgLogElcDensityLocal = 0.0;
    double avgLogElcDensity = 0.0;
    for (int ix = localRgn.getLower(0); ix < localRgn.getUpper(0); ix++)
    {
      nElcIn.setPtr(nElcInPtr, ix);
      Eigen::VectorXd nElcVec(nlocal);

      for (unsigned nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
        nElcVec(nodeIndex) = nElcInPtr[nodeIndex];
      // Calculate nElc at quadrature points
      Eigen::VectorXd nElcAtQuadPoints = interpMatrix*nElcVec;

      for (unsigned componentIndex = 0; componentIndex < nElcAtQuadPoints.rows(); componentIndex++)
      {
        if (nElcAtQuadPoints(componentIndex) != 0.0)
        {
          intElcDensityLocal += gaussWeights[componentIndex]*nElcAtQuadPoints(componentIndex);
          avgLogElcDensityLocal += gaussWeights[componentIndex]*nElcAtQuadPoints(componentIndex)*
            std::log(nElcAtQuadPoints(componentIndex));
        }
      }

      // Initial condition: n_i(z) = n_e(z)
      nIonOut.setPtr(nIonOutPtr, ix);
      for (unsigned nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
        nIonOutPtr[nodeIndex] = nElcInPtr[nodeIndex];

    }
#ifdef HAVE_MPI
    TxCommBase *zComm = grid.getComm();
    zComm->allreduce(1, &intElcDensityLocal, &intElcDensity, TX_SUM);
    zComm->allreduce(1, &avgLogElcDensityLocal, &avgLogElcDensity, TX_SUM);
#else
    intElcDensity    = intElcDensityLocal;
    avgLogElcDensity = avgLogElcDensityLocal;
#endif
    // Divide by n_e
    avgLogElcDensity = avgLogElcDensity/intElcDensity;

    // Loop parameters
    int maxIter = 10000;
    int iter = 0;
    double tol = 1e-14;
    double relError;

    do {
      // Calculate <phi(z)*n_i(z)>
      double ionDensityWeightedPhiLocal = 0.0;
      double ionDensityWeightedPhi = 0.0;
      for (int ix = localRgn.getLower(0); ix < localRgn.getUpper(0); ix++)
      {
        nElcIn.setPtr(nElcInPtr, ix);
        nIonOut.setPtr(nIonOutPtr, ix);

        Eigen::VectorXd nElcVec(nlocal);
        Eigen::VectorXd nIonVec(nlocal);

        for (unsigned nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
        {
          nElcVec(nodeIndex) = nElcInPtr[nodeIndex];
          nIonVec(nodeIndex) = nIonOutPtr[nodeIndex];
        }
       
        Eigen::VectorXd nIonAtQuadPoints = interpMatrix*nIonVec;
        Eigen::VectorXd nElcAtQuadPoints = interpMatrix*nElcVec;
        
        for (unsigned componentIndex = 0; componentIndex < nElcAtQuadPoints.rows(); componentIndex++)
        {
          if (nElcAtQuadPoints(componentIndex) != 0.0)
          {
            double phiAtNode = std::log(nElcAtQuadPoints(componentIndex)) - avgLogElcDensity;
            ionDensityWeightedPhiLocal += gaussWeights[componentIndex]*nIonAtQuadPoints(componentIndex)*phiAtNode;
	  }
        }
      }
#ifdef HAVE_MPI
      zComm->allreduce(1, &ionDensityWeightedPhiLocal, &ionDensityWeightedPhi, TX_SUM);
#else
      ionDensityWeightedPhi = ionDensityWeightedPhiLocal;
#endif
      ionDensityWeightedPhi = ionDensityWeightedPhi/intElcDensity;

      // Loop over all cells to calculate intermediate ion density
      for (int ix = localRgn.getLower(0); ix < localRgn.getUpper(0); ix++)
      {
        nElcIn.setPtr(nElcInPtr, ix);
        nIonOut.setPtr(nIonOutPtr, ix);

        for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
        {
          if (nElcInPtr[nodeIndex] != 0.0)
          {
            double phiAtNode = std::log(nElcInPtr[nodeIndex]) - avgLogElcDensity;
            nIonOutPtr[nodeIndex] = nElcInPtr[nodeIndex]/(1 - 
                kPerpTimesRho*kPerpTimesRho*(phiAtNode - ionDensityWeightedPhi));
          }
        }
      }

      // Calculate total ion number (for <n_i>)
      double intIonDensityLocal = 0.0;
      double intIonDensity = 0.0;
      for (int ix = localRgn.getLower(0); ix < localRgn.getUpper(0); ix++)
      {
        nIonOut.setPtr(nIonOutPtr, ix);

        Eigen::VectorXd nIonVec(nlocal);

        for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
          nIonVec(nodeIndex) = nIonOutPtr[nodeIndex];

        Eigen::VectorXd nIonAtQuadPoints = interpMatrix*nIonVec;

        for (int componentIndex = 0; componentIndex < nIonAtQuadPoints.rows(); componentIndex++)
          intIonDensityLocal += gaussWeights[componentIndex]*nIonAtQuadPoints(componentIndex);
      }
#ifdef HAVE_MPI
      zComm->allreduce(1, &intIonDensityLocal, &intIonDensity, TX_SUM);
#else
      intIonDensity = intIonDensityLocal;
#endif
      
      // Apply correction factor to nIon and calculate relative errors
      double maxResLocal = 0.0;
      double maxRes = 0.0;
      double maxIonDensityLocal = 0.0;
      double maxIonDensity = 0.0;
      for (int ix = localRgn.getLower(0); ix < localRgn.getUpper(0); ix++)
      {
        nIonOut.setPtr(nIonOutPtr, ix);

        for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
        {
          nIonOutPtr[nodeIndex] = nIonOutPtr[nodeIndex] - (intIonDensity - intElcDensity)/domainLength;

          if (fabs(nIonOutPtr[nodeIndex]) > maxIonDensityLocal)
            maxIonDensityLocal = fabs(nIonOutPtr[nodeIndex]);
          // Calculate and possibly store error
          double residual = fabs(nIonOutPtr[nodeIndex] - nIonPrev((ix-localRgn.getLower(0))*nlocal + nodeIndex));
          if (residual > maxResLocal)
            maxResLocal = residual;
          // Store results for comparison
          nIonPrev((ix-localRgn.getLower(0))*nlocal + nodeIndex) = nIonOutPtr[nodeIndex];
        }
      }
#ifdef HAVE_MPI
      // find max residual across all processes.
      zComm->allreduce(1, &maxIonDensityLocal, &maxIonDensity, TX_MAX);
      // find max residual across all processes.
      zComm->allreduce(1, &maxResLocal, &maxRes, TX_MAX);
#else
      maxIonDensity = maxIonDensityLocal;
      maxRes = maxResLocal;
#endif
      relError = maxRes/maxIonDensity;
      std::cout << "iter " << iter << " relerr = " << relError << std::endl;
      iter++;
    } while (iter < maxIter && relError > tol);
 
    return Lucee::UpdaterStatus();
  }

  void
  SOLIonDensityInitialization::declareTypes()
  {
    // takes one input: n_e(x), a 1d DG field
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    // returns one output: n_i(x), a 1d DG field
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  void
  SOLIonDensityInitialization::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}
