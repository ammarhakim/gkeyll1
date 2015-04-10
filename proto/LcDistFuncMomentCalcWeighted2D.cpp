/**
 * @file	LcDistFuncMomentCalcWeighted2D.cpp
 *
 * @brief	Updater to compute 2d moments of a 4d distribution function with an additional weighting function.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDistFuncMomentCalcWeighted2D.h>
#include <LcGlobals.h>
#include <LcLinAlgebra.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <cmath>
#include <vector>

namespace Lucee
{
  const char *DistFuncMomentCalcWeighted2D::id = "DistFuncMomentCalcWeighted2D";

  DistFuncMomentCalcWeighted2D::DistFuncMomentCalcWeighted2D()
    : Lucee::UpdaterIfc()
  {
  }
  
  void
  DistFuncMomentCalcWeighted2D::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    // get hold of 4D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<4> >("basis4d"))
      nodalBasis4d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<4> >("basis4d");
    else
      throw Lucee::Except(
        "DistFuncMomentCalcWeighted2D::readInput: Must specify 4D element to use using 'basis4d'");

    // get hold of 2D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis2d"))
      nodalBasis2d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis2d");
    else
      throw Lucee::Except(
        "DistFuncMomentCalcWeighted2D::readInput: Must specify 2D element to use using 'basis2d'");

    // get moment to compute
    if (tbl.hasNumber("moment"))
    calcMom = (unsigned) tbl.getNumber("moment");
    else
      throw Lucee::Except(
        "DistFuncMomentCalcWeighted2D::readInput: Must specify moment using 'moment'");

    // get moment direction to compute
    if (tbl.hasNumber("momentDirection"))
      momDir = (unsigned) tbl.getNumber("momentDirection");
    else
      momDir = 2;

    if (calcMom > 2)
    {
      Lucee::Except lce("DistFuncMomentCalcWeighted2D::readInput: Only 'moment' < 3 is supported. ");
      lce << "Supplied " << calcMom << " instead";
      throw lce;
    }
  }

  void
  DistFuncMomentCalcWeighted2D::initialize()
  {
    // call base class method
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<4>& grid 
      = this->getGrid<Lucee::StructuredGridBase<4> >();

    // local region to update
    Lucee::Region<4, int> localRgn = grid.getLocalRegion();

    // get number of nodes in 2D and 4D
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

    mom0MatrixVector = std::vector<Eigen::MatrixXd>(nlocal2d);
    mom1MatrixVector = std::vector<Eigen::MatrixXd>(nlocal2d);
    mom2MatrixVector = std::vector<Eigen::MatrixXd>(nlocal2d);

    for (int h = 0; h < nlocal2d; h++)
    {
      mom0MatrixVector[h] = Eigen::MatrixXd::Zero(nlocal2d, nlocal4d);
      mom1MatrixVector[h] = Eigen::MatrixXd::Zero(nlocal2d, nlocal4d);
      mom2MatrixVector[h] = Eigen::MatrixXd::Zero(nlocal2d, nlocal4d);

      for (int i = 0; i < nlocal2d; i++)
      {
        for (int j = 0; j < nlocal4d; j++)
        {
          // Compute integral of phi2d_h * phi2d_i * phi4d_j
          double integralResult[3] = {};
          for (int gaussIndex = 0; gaussIndex < volWeights4d.size(); gaussIndex++)
          {
            double baseIntegral = volWeights4d[gaussIndex]*volQuad2d(gaussIndex % nVolQuad2d, h)*
              volQuad2d(gaussIndex % nVolQuad2d, i)*volQuad4d(gaussIndex, j);
            integralResult[0] += baseIntegral;
            // Get coordinate of quadrautre point in direction momDir
            double coord2Val = volCoords4d(gaussIndex, momDir)*grid.getDx(momDir)/2.0;
            integralResult[1] += coord2Val*baseIntegral;
            integralResult[2] += coord2Val*coord2Val*baseIntegral;
          }
          mom0MatrixVector[h](i, j) = integralResult[0];
          mom1MatrixVector[h](i, j) = integralResult[1];
          mom2MatrixVector[h](i, j) = integralResult[2];
        }
      }
    }

    // Get 2D Mass Matrix
    Lucee::Matrix<double> tempMassMatrix2d(nlocal2d, nlocal2d);
    nodalBasis2d->getMassMatrix(tempMassMatrix2d);
    Eigen::MatrixXd massMatrix2d(nlocal2d, nlocal2d);
    copyLuceeToEigen(tempMassMatrix2d, massMatrix2d);

    // Multiply matrices by inverse of mass matrix
    for (int h = 0; h < nlocal2d; h++)
    {
      mom0MatrixVector[h] = massMatrix2d.inverse()*mom0MatrixVector[h];
      mom1MatrixVector[h] = massMatrix2d.inverse()*mom1MatrixVector[h];
      mom2MatrixVector[h] = massMatrix2d.inverse()*mom2MatrixVector[h];
    }
  }

  Lucee::UpdaterStatus
  DistFuncMomentCalcWeighted2D::update(double t)
  {
    
    // get hold of grid
    const Lucee::StructuredGridBase<4>& grid 
      = this->getGrid<Lucee::StructuredGridBase<4> >();

    // get input field (4d)
    const Lucee::Field<4, double>& distF = this->getInp<Lucee::Field<4, double> >(0);
    // get weighting field (2d)
    const Lucee::Field<2, double>& weightF = this->getInp<Lucee::Field<2, double> >(1);
    // get output field (2d)
    Lucee::Field<2, double>& moment = this->getOut<Lucee::Field<2, double> >(0);

    // local region to update (This is the 4D region. The 2D region is
    // assumed to have the same cell layout as the X-direction of the 4D region)
    Lucee::Region<4, int> localRgn = grid.getLocalRegion();
    Lucee::Region<4, int> localExtRgn = distF.getExtRegion();

    // Make sure we don't integrate over velocity space ghost cells
    localExtRgn.setLower(2, localRgn.getLower(2));
    localExtRgn.setUpper(2, localRgn.getUpper(2));
    localExtRgn.setLower(3, localRgn.getLower(3));
    localExtRgn.setUpper(3, localRgn.getUpper(3));

    // clear out contents of output field
    moment = 0.0;

    // iterators into fields
    Lucee::ConstFieldPtr<double> distFPtr = distF.createConstPtr();
    Lucee::ConstFieldPtr<double> weightFPtr = weightF.createConstPtr();
    Lucee::FieldPtr<double> momentPtr = moment.createPtr();

    int idx[4];
    double xc[4];
    Lucee::RowMajorSequencer<4> seq(localExtRgn);
    unsigned nlocal2d = nodalBasis2d->getNumNodes();
    unsigned nlocal4d = nodalBasis4d->getNumNodes();

    while(seq.step())
    {
      seq.fillWithIndex(idx);

      grid.setIndex(idx);
      grid.getCentroid(xc);

      moment.setPtr(momentPtr, idx[0], idx[1]);
      distF.setPtr(distFPtr, idx);
      weightF.setPtr(weightFPtr, idx[0], idx[1]);

      Eigen::VectorXd distfVec(nlocal4d);
      for (int i = 0; i < nlocal4d; i++)
        distfVec(i) = distFPtr[i];

      // Loop over each component of the weighting function
      for (int h = 0; h < nlocal2d; h++)
      {
        Eigen::VectorXd resultVector(nlocal2d);
        // Calculate contribution to moment
        if (calcMom == 0)
          resultVector = weightFPtr[h]*mom0MatrixVector[h]*distfVec;
        else if (calcMom == 1)
          resultVector = weightFPtr[h]*mom1MatrixVector[h]*distfVec + 
            xc[momDir]*weightFPtr[h]*mom0MatrixVector[h]*distfVec;
        else if (calcMom == 2)
          resultVector = weightFPtr[h]*mom2MatrixVector[h]*distfVec + 
            2*xc[momDir]*weightFPtr[h]*mom1MatrixVector[h]*distfVec +
            xc[momDir]*xc[momDir]*weightFPtr[h]*mom0MatrixVector[h]*distfVec;
        // Accumulate contribution to moment from this cell
        for (int i = 0; i < nlocal2d; i++)
          momentPtr[i] = momentPtr[i] + resultVector(i);
      }
    }

    return Lucee::UpdaterStatus();
  }

  void
  DistFuncMomentCalcWeighted2D::declareTypes()
  {
    // distribution function (4d)
    this->appendInpVarType(typeid(Lucee::Field<4, double>));
    // weighting function (2d)
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    // output field (2d)
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }

  void
  DistFuncMomentCalcWeighted2D::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}
