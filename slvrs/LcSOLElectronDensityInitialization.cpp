/**
 * @file	LcSOLElectronDensityInitialization.cpp
 *
 * @brief	Updater to compute phi using a fixed value of k_perp*rho_s
 * NOTE: CLEAN UP CODE
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLElectronDensityInitialization.h>
#include <LcMathPhysConstants.h>

// gsl includes
#ifdef HAVE_GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#endif

namespace Lucee
{
// set id for module system
  const char *SOLElectronDensityInitialization::id = "SOLElectronDensityInitialization";

  SOLElectronDensityInitialization::SOLElectronDensityInitialization()
    : UpdaterIfc()
  {
  }

  void 
  SOLElectronDensityInitialization::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis");
    else
      throw Lucee::Except("SOLElectronDensityInitialization::readInput: Must specify element to use using 'basis'");

    if (tbl.hasNumber("kPerpTimesRho"))
      kPerpTimesRho = tbl.getNumber("kPerpTimesRho");
    else
      throw Lucee::Except("SOLElectronDensityInitialization::readInput: Must specify kPerpTimesRho");
  }

  void 
  SOLElectronDensityInitialization::initialize()
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
  SOLElectronDensityInitialization::update(double t)
  {
    const Lucee::StructuredGridBase<1>& grid
      = this->getGrid<Lucee::StructuredGridBase<1> >();

    // Electron and ion densities
    const Lucee::Field<1, double>& nIonIn = this->getInp<Lucee::Field<1, double> >(0);
    Lucee::Field<1, double>& nElcOut = this->getOut<Lucee::Field<1, double> >(0);

    int nlocal = nodalBasis->getNumNodes();

    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion();

    Lucee::ConstFieldPtr<double> nIonInPtr = nIonIn.createConstPtr();
    Lucee::FieldPtr<double> nElcOutPtr = nElcOut.createPtr();

    nElcOutPtr = 0.0;

    // Calculate mean ion density
    double avgIonDensity = 0.0;
    for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
    {
      nIonIn.setPtr(nIonInPtr, ix);
      Eigen::VectorXd nIonVec(nlocal);

      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
        nIonVec(componentIndex) = nIonInPtr[componentIndex];
      // Calculate nIon at quadrature points
      Eigen::VectorXd nIonAtQuadPoints = interpMatrix*nIonVec;

      for (int componentIndex = 0; componentIndex < nIonAtQuadPoints.rows(); componentIndex++)
        avgIonDensity += gaussWeights[componentIndex]*nIonAtQuadPoints(componentIndex);
    }
    // Divide by length of domain
    avgIonDensity = avgIonDensity/(grid.getDx(0)*(globalRgn.getUpper(0)-globalRgn.getLower(0)));
    
    int maxIter = 1000;

    // Solve for z by iteration (store temporary results in output structure)
    for(int iter = 0; iter < maxIter; iter++)
    {
      // Calculate <exp(z)>
      double avgExpZ = 0.0;
      for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
      {
        nElcOut.setPtr(nElcOutPtr, ix);
        Eigen::VectorXd expZVec(nlocal);

        for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
          expZVec(componentIndex) = std::exp(nElcOutPtr[componentIndex]);
       
        Eigen::VectorXd expZAtQuadPoints = interpMatrix*expZVec;
        
        for (int componentIndex = 0; componentIndex < expZAtQuadPoints.rows(); componentIndex++)
          avgExpZ += gaussWeights[componentIndex]*expZAtQuadPoints(componentIndex);
      }
      // Divide by length of domain
      avgExpZ = avgExpZ/(grid.getDx(0)*(globalRgn.getUpper(0)-globalRgn.getLower(0)));

      // Loop over all cells
      for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
      {
        // Set inputs
        nIonIn.setPtr(nIonInPtr, ix);
        // Set outputs
        nElcOut.setPtr(nElcOutPtr, ix);

        for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
        {
          double logArg = nIonInPtr[nodeIndex]/avgIonDensity*
             (1-kPerpTimesRho*kPerpTimesRho*nElcOutPtr[nodeIndex])*
             avgExpZ;
          nElcOutPtr[nodeIndex] = std::log(logArg);
        }
      }
    }

    // Calculate <exp(z)>
    double avgExpZ = 0.0;
    for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
    {
      nElcOut.setPtr(nElcOutPtr, ix);
      Eigen::VectorXd expZVec(nlocal);

      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
        expZVec(componentIndex) = std::exp(nElcOutPtr[componentIndex]);
     
      Eigen::VectorXd expZAtQuadPoints = interpMatrix*expZVec;
      
      for (int componentIndex = 0; componentIndex < expZAtQuadPoints.rows(); componentIndex++)
        avgExpZ += gaussWeights[componentIndex]*expZAtQuadPoints(componentIndex);
    }
    // Divide by length of domain
    avgExpZ = avgExpZ/(grid.getDx(0)*(globalRgn.getUpper(0)-globalRgn.getLower(0)));

    // Final loop: take z and compute n_e = <n_i>*exp(z)
    for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
    {
      nElcOut.setPtr(nElcOutPtr, ix);

      for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
        nElcOutPtr[nodeIndex] = avgIonDensity*std::exp(nElcOutPtr[nodeIndex])/avgExpZ;
    }

/*
    // Testing stuff
    const gsl_root_fdfsolver_type *T;
    gsl_root_fdfsolver *s;
    gsl_function_fdf FDF;

    T = gsl_root_fdfsolver_newton;
    s = gsl_root_fdfsolver_alloc (T);

    FDF.f = &Lucee::SOLElectronDensityInitialization::computeF_Wrapper;
    FDF.df = &Lucee::SOLElectronDensityInitialization::computeFPrime_Wrapper;
    FDF.fdf = &Lucee::SOLElectronDensityInitialization::computeFAndFPrime_Wrapper;

    // Loop over all cells
    for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
    {
      // Set inputs
      nIonIn.setPtr(nIonInPtr, ix);
      // Set outputs
      nElcOut.setPtr(nElcOutPtr, ix);

      // Loop over nodes
      for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
      {
        struct fParams params = {meanIonDensity, nIonInPtr[nodeIndex]};
        FDF.params = &params;
        double x0;
        int status;
        double x = (nIonInPtr[nodeIndex]-meanIonDensity)/(meanIonDensity
            + nIonInPtr[nodeIndex]*kPerpTimesRho*kPerpTimesRho);
        double xInit = (nIonInPtr[nodeIndex]-meanIonDensity)/(meanIonDensity
            + nIonInPtr[nodeIndex]*kPerpTimesRho*kPerpTimesRho);
        
        gsl_root_fdfsolver_set (s, &FDF, x);

        int iter = 0;
        int max_iter = 1000;

        do
        {
          iter++;
          status = gsl_root_fdfsolver_iterate (s);
          x0 = x;
          x = gsl_root_fdfsolver_root (s);
          status = gsl_root_test_delta (x, x0, 0, 1e-14);
        }
        while (status == GSL_CONTINUE && iter < max_iter);

        nElcOutPtr[nodeIndex] = meanIonDensity*exp(x);
      }
    }
    
    gsl_root_fdfsolver_free (s);*/

    return Lucee::UpdaterStatus();
  }

  void
  SOLElectronDensityInitialization::declareTypes()
  {
    // takes one input: n_i(x), a 1d DG field
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // returns one output: n_e(x), a 1d DG field
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }

  void
  SOLElectronDensityInitialization::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }

  double
  SOLElectronDensityInitialization::computeF(double x, void *params)
  {
    struct fParams *p = (struct fParams *) params;

    double ionDensityAtPoint = p->ionDensityAtPoint;
    double averageIonDensity = p->averageIonDensity;
    
    return ionDensityAtPoint - averageIonDensity*exp(x) - ionDensityAtPoint*kPerpTimesRho*kPerpTimesRho*x;
  }

  double
  SOLElectronDensityInitialization::computeF_Wrapper(double x, void *params)
  {
    return static_cast<SOLElectronDensityInitialization*>(params)->computeF(x,params);
  }

  double
  SOLElectronDensityInitialization::computeFPrime(double x, void *params)
  {
    struct fParams *p = (struct fParams *) params;

    double ionDensityAtPoint = p->ionDensityAtPoint;
    double averageIonDensity = p->averageIonDensity;
    
    return -averageIonDensity*exp(x) - ionDensityAtPoint*kPerpTimesRho*kPerpTimesRho;
  }

  double
  SOLElectronDensityInitialization::computeFPrime_Wrapper(double x, void *params)
  {
    return static_cast<SOLElectronDensityInitialization*>(params)->computeFPrime(x,params);
  }

  void
  SOLElectronDensityInitialization::computeFAndFPrime(double x, void *params, double *y, double *dy)
  {
    struct fParams *p = (struct fParams *) params;

    double ionDensityAtPoint = p->ionDensityAtPoint;
    double averageIonDensity = p->averageIonDensity;

    *y = ionDensityAtPoint - averageIonDensity*exp(x) - ionDensityAtPoint*kPerpTimesRho*kPerpTimesRho*x;
    *dy = -averageIonDensity*exp(x) - ionDensityAtPoint*kPerpTimesRho*kPerpTimesRho;
  }

  void
  SOLElectronDensityInitialization::computeFAndFPrime_Wrapper(double x, void *params, double *y, double *dy)
  {
    static_cast<SOLElectronDensityInitialization*>(params)->computeFAndFPrime(x,params,y,dy);
  }
}
