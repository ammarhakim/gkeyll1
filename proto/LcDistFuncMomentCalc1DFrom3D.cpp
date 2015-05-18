/**
 * @file	LcDistFuncMomentCalc1DFrom3D.cpp
 *
 * @brief	Updater to compute 1d moments of a 3d distribution function
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDistFuncMomentCalc1DFrom3D.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  const char *DistFuncMomentCalc1DFrom3D::id = "DistFuncMomentCalc1DFrom3D";

  DistFuncMomentCalc1DFrom3D::DistFuncMomentCalc1DFrom3D()
    : Lucee::UpdaterIfc()
  {
  }
  
  void
  DistFuncMomentCalc1DFrom3D::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    // get hold of 3D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except(
        "DistFuncMomentCalc1DFrom3D::readInput: Must specify 3D element to use using 'basis3d'");

    // get hold of 1D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis1d"))
      nodalBasis1d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis1d");
    else
      throw Lucee::Except(
        "DistFuncMomentCalc1DFrom3D::readInput: Must specify 1D element to use using 'basis1d'");

    // get moment to compute
    if (tbl.hasNumber("moment"))
    calcMom = (unsigned) tbl.getNumber("moment");
    else
      throw Lucee::Except(
        "DistFuncMomentCalc1DFrom3D::readInput: Must specify moment using 'moment'");

    // get direction to compute moment in (matters if moment > 1)
    if (tbl.hasNumber("momentDirection"))
      momDir = (unsigned) tbl.getNumber("momentDirection");
    else
    {
      if (calcMom > 0)
      {
      throw Lucee::Except(
        "DistFuncMomentCalc1DFrom3D::readInput: Must specify moment direction using 'momentDirection'");
      }
      else momDir = 1;
    }

    if (calcMom > 2)
    {
      Lucee::Except lce("DistFuncMomentCalc1DFrom3D::readInput: Only 'moment' 0, 1, or 2 is supported. ");
      lce << "Supplied " << calcMom << " instead";
      throw lce;
    }
  }

  void
  DistFuncMomentCalc1DFrom3D::initialize()
  {
    // call base class method
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    // local region to update
    Lucee::Region<3, int> localRgn = grid.getLocalRegion();

    // get number of nodes in 1D and 3D
    unsigned nlocal1d = nodalBasis1d->getNumNodes();
    unsigned nlocal3d = nodalBasis3d->getNumNodes();

    // get volume interpolation matrices for 1d element
    int nVolQuad1d = nodalBasis1d->getNumGaussNodes();
    std::vector<double> volWeights1d(nVolQuad1d);
    Lucee::Matrix<double> tempVolQuad1d(nVolQuad1d, nlocal1d);
    Lucee::Matrix<double> tempVolCoords1d(nVolQuad1d, NC1);

    nodalBasis1d->getGaussQuadData(tempVolQuad1d, tempVolCoords1d, volWeights1d);

    Eigen::MatrixXd volQuad1d(nVolQuad1d, nlocal1d);
    copyLuceeToEigen(tempVolQuad1d, volQuad1d);
    // TESTING STUFF
    Eigen::MatrixXd volCoords1d(nVolQuad1d, NC1);
    copyLuceeToEigen(tempVolCoords1d, volCoords1d);

    // get volume interpolation matrices for 3d element
    int nVolQuad3d = nodalBasis3d->getNumGaussNodes();
    std::vector<double> volWeights3d(nVolQuad3d);
    Lucee::Matrix<double> tempVolQuad3d(nVolQuad3d, nlocal3d);
    Lucee::Matrix<double> tempVolCoords3d(nVolQuad3d, NC3);

    nodalBasis3d->getGaussQuadData(tempVolQuad3d, tempVolCoords3d, volWeights3d);

    Eigen::MatrixXd volQuad3d(nVolQuad3d, nlocal3d);
    copyLuceeToEigen(tempVolQuad3d, volQuad3d);
    // TESTING STUFF
    Eigen::MatrixXd volCoords3d(nVolQuad3d, NC3);
    copyLuceeToEigen(tempVolCoords3d, volCoords3d);

    mom0Matrix = Eigen::MatrixXd::Zero(nlocal1d, nlocal3d);
    mom1Matrix = Eigen::MatrixXd::Zero(nlocal1d, nlocal3d);
    mom2Matrix = Eigen::MatrixXd::Zero(nlocal1d, nlocal3d);

    for (int i = 0; i < nlocal1d; i++)
    {
      for (int j = 0; j < nlocal3d; j++)
      {
        // Compute integral of phi2d_i * phi4d_j
        double integralResult[3] = {};
        for (int gaussIndex = 0; gaussIndex < volWeights3d.size(); gaussIndex++)
        {
          double baseIntegral = volWeights3d[gaussIndex]*volQuad1d(gaussIndex % nVolQuad1d, i)*
            volQuad3d(gaussIndex, j);
          integralResult[0] += baseIntegral;
          // Get coordinate of quadrature point in direction momDir
          double coord2Val = volCoords3d(gaussIndex, momDir)*grid.getDx(momDir)/2.0;
          integralResult[1] += coord2Val*baseIntegral;
          integralResult[2] += coord2Val*coord2Val*baseIntegral;
        }
        mom0Matrix(i, j) = integralResult[0];
        mom1Matrix(i, j) = integralResult[1];
        mom2Matrix(i, j) = integralResult[2];
      }
    }

    // Get 1D Mass Matrix
    Lucee::Matrix<double> tempMassMatrix1d(nlocal1d, nlocal1d);
    nodalBasis1d->getMassMatrix(tempMassMatrix1d);
    Eigen::MatrixXd massMatrix1d(nlocal1d, nlocal1d);
    copyLuceeToEigen(tempMassMatrix1d, massMatrix1d);

    // Multiply matrices by inverse of mass matrix
    mom0Matrix = massMatrix1d.inverse()*mom0Matrix;
    mom1Matrix = massMatrix1d.inverse()*mom1Matrix;
    mom2Matrix = massMatrix1d.inverse()*mom2Matrix;
  }

  Lucee::UpdaterStatus
  DistFuncMomentCalc1DFrom3D::update(double t)
  {
    
    // get hold of grid
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    // get input field (3d)
    const Lucee::Field<3, double>& distF = this->getInp<Lucee::Field<3, double> >(0);
    // get output field (1d)
    Lucee::Field<1, double>& moment = this->getOut<Lucee::Field<1, double> >(0);

    // local region to update (This is the 3D region. The 1D region is
    // assumed to have the same cell layout as the X-direction of the 3D region)
    Lucee::Region<3, int> localRgn = grid.getLocalRegion();

    // clear out contents of output field
    moment = 0.0;

    // iterators into fields
    Lucee::ConstFieldPtr<double> distFPtr = distF.createConstPtr();
    Lucee::FieldPtr<double> momentPtr = moment.createPtr();

    int idx[3];
    double xc[3];
    Lucee::RowMajorSequencer<3> seq(localRgn);
    unsigned nlocal1d = nodalBasis1d->getNumNodes();
    unsigned nlocal3d = nodalBasis3d->getNumNodes();

    while(seq.step())
    {
      seq.fillWithIndex(idx);

      grid.setIndex(idx);
      grid.getCentroid(xc);

      moment.setPtr(momentPtr, idx[0]);
      distF.setPtr(distFPtr, idx);

      Eigen::VectorXd distfVec(nlocal3d);
      for (int i = 0; i < nlocal3d; i++)
        distfVec(i) = distFPtr[i];

      // Accumulate contribution to moment from this cell
      Eigen::VectorXd resultVector(nlocal1d);

      if (calcMom == 0)
        resultVector = mom0Matrix*distfVec;
      else if (calcMom == 1)
        resultVector = mom1Matrix*distfVec + xc[momDir]*mom0Matrix*distfVec;
      else if (calcMom == 2)
        resultVector = mom2Matrix*distfVec + 2*xc[momDir]*mom1Matrix*distfVec +
          xc[momDir]*xc[momDir]*mom0Matrix*distfVec;

      for (int i = 0; i < nlocal1d; i++)
        momentPtr[i] = momentPtr[i] + resultVector(i);
    }

    return Lucee::UpdaterStatus();
  }

  void
  DistFuncMomentCalc1DFrom3D::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }

  void
  DistFuncMomentCalc1DFrom3D::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}
