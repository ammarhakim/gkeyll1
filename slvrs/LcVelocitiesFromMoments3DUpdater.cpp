/**
 * @file	LcVelocitiesFromMoments3DUpdater.cpp
 *
 * @brief	Updater to compute parallel drift velocity and thermal velocity squared
 * given certain moments of the distribution function. Used for 3D SOL problem.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAlignedRectCoordSys.h>
#include <LcField.h>
#include <LcLinAlgebra.h>
#include <LcMathLib.h>
#include <LcVelocitiesFromMoments3DUpdater.h>
#include <LcStructuredGridBase.h>

// math include
#include <cmath>

namespace Lucee
{
// set id for module system
  const char *VelocitiesFromMoments3DUpdater::id = "VelocitiesFromMoments3DUpdater";

  VelocitiesFromMoments3DUpdater::VelocitiesFromMoments3DUpdater()
    : UpdaterIfc()
  {
  }

  void 
  VelocitiesFromMoments3DUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis");
    else
      throw Lucee::Except("VelocitiesFromMoments3DUpdater::readInput: Must specify element to use using 'basis1d'");
  }

  void 
  VelocitiesFromMoments3DUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();
    // local region to update
    Lucee::Region<1, int> localRgn = grid.getLocalRegion();

    Lucee::RowMajorSequencer<1> seq(localRgn);
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

    gaussWeights = std::vector<double>(numQuadNodes);

    // Allocate Eigen matrices
    interpMatrix = Eigen::MatrixXd(numQuadNodes, nlocal);
    gaussOrdinates = Eigen::MatrixXd(numQuadNodes, 3);

    // Get the interpolation matrix for the volume quadrature points
    nodalBasis->getGaussQuadData(interpMatrixLucee, gaussOrdinatesLucee, gaussWeights);

    copyLuceeToEigen(massMatrixLucee, massMatrix);
    copyLuceeToEigen(interpMatrixLucee, interpMatrix);
    copyLuceeToEigen(gaussOrdinatesLucee, gaussOrdinates);

    // Compute and store inverse of mass matrix
    massMatrixInv = massMatrix.inverse();

    // Compute and store transpose of interpolation matrix
    interpMatrixTranspose = massMatrixInv*interpMatrix.transpose();
  }

  Lucee::UpdaterStatus 
  VelocitiesFromMoments3DUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<1>& grid
      = this->getGrid<Lucee::StructuredGridBase<1> >();

    // Particle density n(x) = <1>
    const Lucee::Field<1, double>& mom0In = this->getInp<Lucee::Field<1, double> >(0);
    // Drift velocity <v_para>
    const Lucee::Field<1, double>& mom1ParaIn = this->getInp<Lucee::Field<1, double> >(1);
    // Perpendicular velocity <mu>
    const Lucee::Field<1, double>& mom1MuIn = this->getInp<Lucee::Field<1, double> >(2);
    // Second velocity moment <v_para^2>(x)
    const Lucee::Field<1, double>& mom2ParaIn = this->getInp<Lucee::Field<1, double> >(3);
    // Drift velocity u_para(x)
    Lucee::Field<1, double>& uOut = this->getOut<Lucee::Field<1, double> >(0);
    // Thermal velocity squared vt(x)^2
    Lucee::Field<1, double>& vtSqOut = this->getOut<Lucee::Field<1, double> >(1);

    int nlocal = nodalBasis->getNumNodes();

    Lucee::Region<1, int> localRgn = grid.getLocalRegion();

    Lucee::ConstFieldPtr<double> mom0Ptr = mom0In.createConstPtr();
    Lucee::ConstFieldPtr<double> mom1ParaPtr = mom1ParaIn.createConstPtr();
    Lucee::ConstFieldPtr<double> mom1MuPtr = mom1MuIn.createConstPtr();
    Lucee::ConstFieldPtr<double> mom2ParaPtr = mom2ParaIn.createConstPtr();
    Lucee::FieldPtr<double> uPtr = uOut.createPtr();
    Lucee::FieldPtr<double> vtSqPtr = vtSqOut.createPtr();

    uPtr = 0.0;
    vtSqPtr = 0.0;

    for (int ix = localRgn.getLower(0); ix < localRgn.getUpper(0); ix++)
    {
      // Set inputs
      mom0In.setPtr(mom0Ptr, ix);
      mom1ParaIn.setPtr(mom1ParaPtr, ix);
      mom1MuIn.setPtr(mom1MuPtr, ix);
      mom2ParaIn.setPtr(mom2ParaPtr, ix);
      // Set outputs
      uOut.setPtr(uPtr, ix);
      vtSqOut.setPtr(vtSqPtr, ix);


      Eigen::VectorXd uVec(nlocal);
      // Compute u(x) naively by dividing weights of n*u by n(x)
      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
      {
        if (mom0Ptr[componentIndex] == 0.0)
        {
          uPtr[componentIndex] = 0.0;
          vtSqPtr[componentIndex] = 0.0;
        }
        else 
        {
          uPtr[componentIndex] = mom1ParaPtr[componentIndex]/mom0Ptr[componentIndex];
          // Fill in first part of vt(x)^2
          vtSqPtr[componentIndex] = (mom1MuPtr[componentIndex] + mom2ParaPtr[componentIndex])/
            (3.0*mom0Ptr[componentIndex]);
        }
        uVec(componentIndex) = uPtr[componentIndex];
      }

      // Compute projection of u(x)^2/3.0
      Eigen::VectorXd uAtQuadPoints = interpMatrix*uVec;
      Eigen::VectorXd uAtQuadPointsSq = Eigen::VectorXd(interpMatrix.rows());

      // Compute u^2*weight/3.0 at each quadrature point
      for (int componentIndex = 0; componentIndex < uAtQuadPointsSq.rows(); componentIndex++)
        uAtQuadPointsSq(componentIndex) = gaussWeights[componentIndex]*uAtQuadPoints(componentIndex)*uAtQuadPoints(componentIndex)/3.0;

      Eigen::VectorXd uSqWeights = interpMatrixTranspose*uAtQuadPoints;

      // Fill in second part of vt(x)^2
      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
        vtSqPtr[componentIndex] -= uSqWeights(componentIndex);
    }

    return Lucee::UpdaterStatus();
  }

  void
  VelocitiesFromMoments3DUpdater::declareTypes()
  {
    // takes four inputs (moment0=<1>=n, moment1=<v_para>, 
    // moment1=<mu>, moment2=<v_para^2>) 
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // returns two outputs (u_para(x), vt^2(x))
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }

  void
  VelocitiesFromMoments3DUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}
