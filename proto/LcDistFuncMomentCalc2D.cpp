/**
 * @file	LcDistFuncMomentCalc2D.h
 *
 * @brief	Updater to compute moments of distribution function.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDistFuncMomentCalc2D.h>
#include <LcGlobals.h>
#include <LcLinAlgebra.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <cmath>
#include <vector>

namespace Lucee
{
  const char *DistFuncMomentCalc2D::id = "DistFuncMomentCalc2D";

  DistFuncMomentCalc2D::DistFuncMomentCalc2D()
    : Lucee::UpdaterIfc()
  {
  }
  
  void
  DistFuncMomentCalc2D::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    // get hold of 4D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<4> >("basis4d"))
      nodalBasis4d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<4> >("basis4d");
    else
      throw Lucee::Except(
        "DistFuncMomentCalc2D::readInput: Must specify 4D element to use using 'basis4d'");

    // get hold of 2D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis2d"))
      nodalBasis2d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis2d");
    else
      throw Lucee::Except(
        "DistFuncMomentCalc2D::readInput: Must specify 2D element to use using 'basis1d'");

    // get moment to compute
    if (tbl.hasNumber("moment"))
    calcMom = (unsigned) tbl.getNumber("moment");
    else
      throw Lucee::Except(
        "DistFuncMomentCalc2D::readInput: Must specify moment using 'moment'");

    if (calcMom > 3)
    {
      Lucee::Except lce("DistFuncMomentCalc2D::readInput: Only 'moment' 0, 1 2, or 3 is supported. ");
      lce << "Supplied " << calcMom << " instead";
      throw lce;
    }
  }

  void
  DistFuncMomentCalc2D::initialize()
  {
    // call base class method
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<4>& grid 
      = this->getGrid<Lucee::StructuredGridBase<4> >();

    // local region to update
    Lucee::Region<4, int> localRgn = grid.getLocalRegion();

    // get number of nodes in 1D and 2D
    unsigned nlocal2d = nodalBasis2d->getNumNodes();
    unsigned nlocal4d = nodalBasis4d->getNumNodes();

    // get volume interpolation matrices for 2d element
    int nVolQuad2d = nodalBasis2d->getNumGaussNodes();
    std::vector<double> volWeights2d(nVolQuad2d);
    Lucee::Matrix<double> tempVolQuad2d(nVolQuad2d, nlocal2d);
    Lucee::Matrix<double> tempVolCoords2d(nVolQuad2d, 3);

    nodalBasis2d->getGaussQuadData(tempVolQuad2d, tempVolCoords2d, volWeights2d);

    Eigen::MatrixXd volQuad2d(nVolQuad2d, nlocal2d);
    copyLuceeToEigen(tempVolQuad2d, volQuad2d);
    // TESTING STUFF
    Eigen::MatrixXd volCoords2d(nVolQuad2d, 3);
    copyLuceeToEigen(tempVolCoords2d, volCoords2d);

    // get volume interpolation matrices for 4d element
    int nVolQuad4d = nodalBasis4d->getNumGaussNodes();
    std::vector<double> volWeights4d(nVolQuad4d);
    Lucee::Matrix<double> tempVolQuad4d(nVolQuad4d, nlocal4d);
    Lucee::Matrix<double> tempVolCoords4d(nVolQuad4d, 4);

    nodalBasis4d->getGaussQuadData(tempVolQuad4d, tempVolCoords4d, volWeights4d);

    Eigen::MatrixXd volQuad4d(nVolQuad4d, nlocal4d);
    copyLuceeToEigen(tempVolQuad4d, volQuad4d);
    // TESTING STUFF
    Eigen::MatrixXd volCoords4d(nVolQuad4d, 4);
    copyLuceeToEigen(tempVolCoords4d, volCoords4d);

    mom0Matrix = Eigen::MatrixXd(nlocal2d, nlocal4d);
    // NOTE: be careful when computing these matrices. need to convert
    // v to physical space v
    mom1Matrix = Eigen::MatrixXd(nlocal2d, nlocal4d);
    mom2Matrix = Eigen::MatrixXd(nlocal2d, nlocal4d);

    double dV1 = grid.getDx(2);
    double dV2 = grid.getDx(3);

    for (int i = 0; i < nlocal2d; i++)
    {
      for (int j = 0; j < nlocal4d; j++)
      {
        // Compute integral of phi2d_i * phi4d_j
        double integralResult = 0.0;
        for (int gaussIndex = 0; gaussIndex < volWeights4d.size(); gaussIndex++)
          integralResult += volWeights4d[gaussIndex]*volQuad2d(gaussIndex % nVolQuad2d, i)*volQuad4d(gaussIndex, j);
        mom0Matrix(i, j) = integralResult;
      }
    }

    // Get 2D Mass Matrix
    Lucee::Matrix<double> tempMassMatrix2d(nlocal2d, nlocal2d);
    nodalBasis2d->getMassMatrix(tempMassMatrix2d);
    Eigen::MatrixXd massMatrix2d(nlocal2d, nlocal2d);
    copyLuceeToEigen(tempMassMatrix2d, massMatrix2d);

    // Multiply matrices by inverse of mass matrix
    mom0Matrix = massMatrix2d.inverse()*mom0Matrix;
  }

  Lucee::UpdaterStatus
  DistFuncMomentCalc2D::update(double t)
  {
    
    // get hold of grid
    const Lucee::StructuredGridBase<4>& grid 
      = this->getGrid<Lucee::StructuredGridBase<4> >();

    // get input field (4d)
    const Lucee::Field<4, double>& distF = this->getInp<Lucee::Field<4, double> >(0);
    // get output field (2D)
    Lucee::Field<2, double>& moment = this->getOut<Lucee::Field<2, double> >(0);

  // local region to update (This is the 4D region. The 2D region is
  // assumed to have the same cell layout as the X-direction of the 4D region)
    Lucee::Region<4, int> localRgn = grid.getLocalRegion();

    // clear out contents of output field
    moment = 0.0;

    // iterators into fields
    Lucee::ConstFieldPtr<double> distFPtr = distF.createConstPtr();
    Lucee::FieldPtr<double> momentPtr = moment.createPtr();

    //double dv = grid.getDx(1);
    //double dv2 = 0.5*dv;
    //double dv22 = dv2*dv2;
    //double dv23 = dv2*dv2*dv2;

    int idx[4];
    double xc[4];
    Lucee::RowMajorSequencer<4> seq(localRgn);
    unsigned nlocal2d = nodalBasis2d->getNumNodes();
    unsigned nlocal4d = nodalBasis4d->getNumNodes();

    while(seq.step())
    {
      seq.fillWithIndex(idx);

      grid.setIndex(idx);
      //grid.getCentroid(xc);

      moment.setPtr(momentPtr, idx[0], idx[1]);
      distF.setPtr(distFPtr, idx);

      Eigen::VectorXd distfVec(nlocal4d);
      for (int i = 0; i < nlocal4d; i++)
        distfVec(i) = distFPtr[i];

      // Accumulate contribution to moment from this cell
      Eigen::VectorXd resultVector(nlocal2d);
      if (calcMom == 0)
        resultVector = mom0Matrix*distfVec;

      for (int i = 0; i < nlocal2d; i++)
        momentPtr[i] = momentPtr[i] + resultVector(i);
    }

    return Lucee::UpdaterStatus();
  }

  void
  DistFuncMomentCalc2D::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<4, double>));
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }

  void 
  DistFuncMomentCalc2D::matVec(double m, const Lucee::Matrix<double>& mat,
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
  DistFuncMomentCalc2D::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}
