/**
 * @file	LcDistFuncMomentCalc3D.cpp
 *
 * @brief	Updater to compute 3d moments of a 5d distribution function.
 * Currently only works for 0th moment
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDistFuncMomentCalc3D.h>
#include <LcGlobals.h>
#include <LcLinAlgebra.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <cmath>
#include <vector>

namespace Lucee
{
  const char *DistFuncMomentCalc3D::id = "DistFuncMomentCalc3D";

  DistFuncMomentCalc3D::DistFuncMomentCalc3D()
    : Lucee::UpdaterIfc()
  {
  }
  
  void
  DistFuncMomentCalc3D::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    // get hold of 4D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except(
        "DistFuncMomentCalc3D::readInput: Must specify 5D element to use using 'basis5d'");

    // get hold of 2D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except(
        "DistFuncMomentCalc3D::readInput: Must specify 3D element to use using 'basis1d'");

    // get moment to compute
    if (tbl.hasNumber("moment"))
    calcMom = (unsigned) tbl.getNumber("moment");
    else
      throw Lucee::Except(
        "DistFuncMomentCalc3D::readInput: Must specify moment using 'moment'");

    if (calcMom > 3)
    {
      Lucee::Except lce("DistFuncMomentCalc3D::readInput: Only 'moment' 0, 1 2, or 3 is supported. ");
      lce << "Supplied " << calcMom << " instead";
      throw lce;
    }
  }

  void
  DistFuncMomentCalc3D::initialize()
  {
    // call base class method
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // local region to update
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();

    // get number of nodes in 1D and 2D
    unsigned nlocal3d = nodalBasis3d->getNumNodes();
    unsigned nlocal5d = nodalBasis5d->getNumNodes();

    // get volume interpolation matrices for 3d element
    int nVolQuad3d = nodalBasis3d->getNumGaussNodes();
    std::vector<double> volWeights3d(nVolQuad3d);
    Lucee::Matrix<double> tempVolQuad3d(nVolQuad3d, nlocal3d);
    Lucee::Matrix<double> tempVolCoords3d(nVolQuad3d, 3);

    nodalBasis3d->getGaussQuadData(tempVolQuad3d, tempVolCoords3d, volWeights3d);

    Eigen::MatrixXd volQuad3d(nVolQuad3d, nlocal3d);
    copyLuceeToEigen(tempVolQuad3d, volQuad3d);
    // TESTING STUFF
    Eigen::MatrixXd volCoords3d(nVolQuad3d, 3);
    copyLuceeToEigen(tempVolCoords3d, volCoords3d);

    // get volume interpolation matrices for 5d element
    int nVolQuad5d = nodalBasis5d->getNumGaussNodes();
    std::vector<double> volWeights5d(nVolQuad5d);
    Lucee::Matrix<double> tempVolQuad5d(nVolQuad5d, nlocal5d);
    Lucee::Matrix<double> tempVolCoords5d(nVolQuad5d, 5);

    nodalBasis5d->getGaussQuadData(tempVolQuad5d, tempVolCoords5d, volWeights5d);

    Eigen::MatrixXd volQuad5d(nVolQuad5d, nlocal5d);
    copyLuceeToEigen(tempVolQuad5d, volQuad5d);
    // TESTING STUFF
    Eigen::MatrixXd volCoords5d(nVolQuad5d, 5);
    copyLuceeToEigen(tempVolCoords5d, volCoords5d);

    mom0Matrix = Eigen::MatrixXd(nlocal3d, nlocal5d);
    // NOTE: be careful when computing these matrices. need to convert
    // v to physical space v
    mom1Matrix = Eigen::MatrixXd(nlocal3d, nlocal5d);
    mom2Matrix = Eigen::MatrixXd(nlocal3d, nlocal5d);

    for (int i = 0; i < nlocal3d; i++)
    {
      for (int j = 0; j < nlocal5d; j++)
      {
        // Compute integral of phi3d_i * phi5d_j
        double integralResult = 0.0;
        for (int gaussIndex = 0; gaussIndex < volWeights5d.size(); gaussIndex++)
          integralResult += volWeights5d[gaussIndex]*volQuad3d(gaussIndex % nVolQuad3d, i)*volQuad5d(gaussIndex, j);
        mom0Matrix(i, j) = integralResult;
      }
    }

    // Get 3D Mass Matrix
    Lucee::Matrix<double> tempMassMatrix3d(nlocal3d, nlocal3d);
    nodalBasis3d->getMassMatrix(tempMassMatrix3d);
    Eigen::MatrixXd massMatrix3d(nlocal3d, nlocal3d);
    copyLuceeToEigen(tempMassMatrix3d, massMatrix3d);

    // Multiply matrices by inverse of mass matrix
    mom0Matrix = massMatrix3d.inverse()*mom0Matrix;
  }

  Lucee::UpdaterStatus
  DistFuncMomentCalc3D::update(double t)
  {
    
    // get hold of grid
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // get input field (5d)
    const Lucee::Field<5, double>& distF = this->getInp<Lucee::Field<5, double> >(0);
    // get output field (3D)
    Lucee::Field<3, double>& moment = this->getOut<Lucee::Field<3, double> >(0);

    // local region to update (This is the 4D region. The 2D region is
    // assumed to have the same cell layout as the X-direction of the 4D region)
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();

    // clear out contents of output field
    moment = 0.0;

    // iterators into fields
    Lucee::ConstFieldPtr<double> distFPtr = distF.createConstPtr();
    Lucee::FieldPtr<double> momentPtr = moment.createPtr();

    int idx[5];
    double xc[5];
    Lucee::RowMajorSequencer<5> seq(localRgn);
    unsigned nlocal3d = nodalBasis3d->getNumNodes();
    unsigned nlocal5d = nodalBasis5d->getNumNodes();

    while(seq.step())
    {
      seq.fillWithIndex(idx);

      grid.setIndex(idx);
      //grid.getCentroid(xc);

      moment.setPtr(momentPtr, idx[0], idx[1], idx[2]);
      distF.setPtr(distFPtr, idx);

      Eigen::VectorXd distfVec(nlocal5d);
      for (int i = 0; i < nlocal5d; i++)
        distfVec(i) = distFPtr[i];

      // Accumulate contribution to moment from this cell
      Eigen::VectorXd resultVector(nlocal3d);
      if (calcMom == 0)
        resultVector = mom0Matrix*distfVec;

      for (int i = 0; i < nlocal3d; i++)
        momentPtr[i] = momentPtr[i] + resultVector(i);
    }

    return Lucee::UpdaterStatus();
  }

  void
  DistFuncMomentCalc3D::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
  }

  void 
  DistFuncMomentCalc3D::matVec(double m, const Lucee::Matrix<double>& mat,
    const double* vec, double v, double *out)
  {
    double tv;
    unsigned rows = mat.numRows(), cols = mat.numColumns();
    for (unsigned i=0; i<rows; ++i)
    {
      tv = 0.0;
      for (unsigned j=0; j<cols; ++j)
        tv += mat(i,j)*vec[j];
      out[i] = m*tv + v*out[i];
    }
  }

  void
  DistFuncMomentCalc3D::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}
