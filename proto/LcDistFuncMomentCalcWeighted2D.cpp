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

  DistFuncMomentCalcWeighted2D::~DistFuncMomentCalcWeighted2D()
  {
    delete moment;
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
            // Get coordinate of quadrature point in direction momDir
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

    int lower[2], upper[2];
    lower[0] = localRgn.getLower(0); upper[0] = localRgn.getUpper(0);
    lower[1] = localRgn.getLower(1); upper[1] = localRgn.getUpper(1);
    Lucee::Region<2, int> local2D(lower, upper);
    int lg[2] = {1,1}, ug[2] = {1,1}; // one ghost cell layer in each direction
    // allocate space for storing local moment calculation
    moment = new Lucee::Field<2, double>(local2D, nlocal2d, lg, ug);
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
    Lucee::Field<2, double>& momentOut = this->getOut<Lucee::Field<2, double> >(0);

    // create duplicate to store local moments
    //Lucee::Field<2, double> moment = momentOut.duplicate();

    // local region to update (This is the 4D region. The 2D region is
    // assumed to have the same cell layout as the X-direction of the 4D region)
    Lucee::Region<4, int> localRgn = grid.getLocalRegion();
    Lucee::Region<4, int> localExtRgn = distF.getExtRegion();

    // Make sure we integrate over conf. space ghost cells
    localRgn.setLower(0, localExtRgn.getLower(0));
    localRgn.setUpper(0, localExtRgn.getUpper(0));
    localRgn.setLower(1, localExtRgn.getLower(1));
    localRgn.setUpper(1, localExtRgn.getUpper(1));
    
    // clear out contents of output field
    (*moment) = 0.0;

    // iterators into fields
    Lucee::ConstFieldPtr<double> distFPtr = distF.createConstPtr();
    Lucee::ConstFieldPtr<double> weightFPtr = weightF.createConstPtr();
    Lucee::FieldPtr<double> momentPtr = moment->createPtr();

    int idx[4];
    double xc[4];
    Lucee::RowMajorSequencer<4> seq(localRgn);
    //Lucee::RowMajorIndexer<4> idxr(localRgn);
    unsigned nlocal2d = nodalBasis2d->getNumNodes();
    unsigned nlocal4d = nodalBasis4d->getNumNodes();

    int lower2d[] = {localRgn.getLower(0), localRgn.getLower(1)};
    int upper2d[] = {localRgn.getUpper(0), localRgn.getUpper(1)};
    Lucee::Region<2, int> rgn2d(lower2d, upper2d);
    Lucee::RowMajorIndexer<2> idxr(rgn2d);
    Lucee::RowMajorSequencer<2> seq2d(rgn2d);

    int localPositionCells = localRgn.getShape(0)*localRgn.getShape(1);
    std::vector<double> localMoment(localPositionCells*nlocal2d);
    Eigen::VectorXd resultVector(nlocal2d);
    Eigen::VectorXd distfVec(nlocal4d);

    while(seq.step())
    {
      seq.fillWithIndex(idx);
      grid.setIndex(idx);
      grid.getCentroid(xc);

      moment->setPtr(momentPtr, idx[0], idx[1]);
      distF.setPtr(distFPtr, idx);
      weightF.setPtr(weightFPtr, idx[0], idx[1]);

      for (int i = 0; i < nlocal4d; i++)
        distfVec(i) = distFPtr[i];

      // Loop over each component of the weighting function
      for (int h = 0; h < nlocal2d; h++)
      {
        // Calculate contribution to moment
        if (calcMom == 0)
          resultVector.noalias() = weightFPtr[h]*mom0MatrixVector[h]*distfVec;
        else if (calcMom == 1)
          resultVector.noalias() = (weightFPtr[h]*mom1MatrixVector[h] + 
            xc[momDir]*weightFPtr[h]*mom0MatrixVector[h])*distfVec;
        else if (calcMom == 2)
          resultVector.noalias() = (weightFPtr[h]*mom2MatrixVector[h] + 
            2*xc[momDir]*weightFPtr[h]*mom1MatrixVector[h] +
            xc[momDir]*xc[momDir]*weightFPtr[h]*mom0MatrixVector[h])*distfVec;
        // Accumulate contribution to moment from this cell
        for (int i = 0; i < nlocal2d; i++)
          momentPtr[i] = momentPtr[i] + resultVector(i);
      }
    }

    int idx2d[2];
    // Loop over each 'local' position space cell
    while(seq2d.step())
    {
      seq2d.fillWithIndex(idx2d);
      int cellIndex = idxr.getIndex(idx2d);

      moment->setPtr(momentPtr, idx2d);
      // copy data to vector
      for (int i = 0; i < nlocal2d; i++)
        localMoment[cellIndex*nlocal2d+i] = momentPtr[i];
    }

    // Above loop computes moments on local phase-space domain. We need to
    // sum across velocity space to get total moment on configuration
    // space.
    std::vector<double> reducedMoment(localPositionCells*nlocal2d);
    // we need to get moment communicator of field as updater's moment
    // communicator is same as its grid's moment communicator. In this
    // case, grid is phase-space grid, which is not what we want.
    TxCommBase *momComm = momentOut.getMomComm();
    // amount to communicate
    momComm->allreduce(localPositionCells*nlocal2d, localMoment, reducedMoment, TX_SUM);

    seq2d.reset();
    Lucee::FieldPtr<double> momentOutPtr = momentOut.createPtr();
    // Copy reducedMoment to output field
    while(seq2d.step())
    {
      seq2d.fillWithIndex(idx2d);
      int cellIndex = idxr.getIndex(idx2d);
      momentOut.setPtr(momentOutPtr, idx2d);
      // copy data to vector
      for (int i = 0; i < nlocal2d; i++)
        momentOutPtr[i] = reducedMoment[cellIndex*nlocal2d+i];
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
