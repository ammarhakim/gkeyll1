/**
 * @file	LcDistFuncMomentCalcWeighted3D.cpp
 *
 * @brief	Updater to compute 3d moments of a 5d distribution function with an additional weighting function.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDistFuncMomentCalcWeighted3D.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <cmath>
#include <vector>

namespace Lucee
{
  const char *DistFuncMomentCalcWeighted3D::id = "DistFuncMomentCalcWeighted3D";

  DistFuncMomentCalcWeighted3D::DistFuncMomentCalcWeighted3D()
    : Lucee::UpdaterIfc()
  {
  }
  
  void
  DistFuncMomentCalcWeighted3D::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    // get hold of 5D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except(
        "DistFuncMomentCalcWeighted3D::readInput: Must specify 5D element to use using 'basis5d'");

    // get hold of 3D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except(
        "DistFuncMomentCalcWeighted3D::readInput: Must specify 3D element to use using 'basis3d'");

    // get moment to compute
    if (tbl.hasNumber("moment"))
    calcMom = (unsigned) tbl.getNumber("moment");
    else
      throw Lucee::Except(
        "DistFuncMomentCalcWeighted3D::readInput: Must specify moment using 'moment'");

    // get moment direction to compute
    if (tbl.hasNumber("momentDirection"))
      momDir = (unsigned) tbl.getNumber("momentDirection");
    else
      momDir = 2;

    if (calcMom > 2)
    {
      Lucee::Except lce("DistFuncMomentCalcWeighted3D::readInput: Only 'moment' < 3 is supported. ");
      lce << "Supplied " << calcMom << " instead";
      throw lce;
    }
  }

  void
  DistFuncMomentCalcWeighted3D::initialize()
  {
    // call base class method
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // local region to update
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();

    // get number of nodes in 3D and 5D
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

    mom0MatrixVector = std::vector<Eigen::MatrixXd>(nlocal3d);
    mom1MatrixVector = std::vector<Eigen::MatrixXd>(nlocal3d);
    mom2MatrixVector = std::vector<Eigen::MatrixXd>(nlocal3d);

    for (int h = 0; h < nlocal3d; h++)
    {
      mom0MatrixVector[h] = Eigen::MatrixXd::Zero(nlocal3d, nlocal5d);
      mom1MatrixVector[h] = Eigen::MatrixXd::Zero(nlocal3d, nlocal5d);
      mom2MatrixVector[h] = Eigen::MatrixXd::Zero(nlocal3d, nlocal5d);

      for (int i = 0; i < nlocal3d; i++)
      {
        for (int j = 0; j < nlocal5d; j++)
        {
          // Compute integral of phi3d_h * phi3d_i * phi5d_j
          double integralResult[3] = {};
          for (int gaussIndex = 0; gaussIndex < volWeights5d.size(); gaussIndex++)
          {
            double baseIntegral = volWeights5d[gaussIndex]*volQuad3d(gaussIndex % nVolQuad3d, h)*
              volQuad3d(gaussIndex % nVolQuad3d, i)*volQuad5d(gaussIndex, j);
            integralResult[0] += baseIntegral;
            // Get coordinate of quadrature point in direction momDir
            double coord3Val = volCoords5d(gaussIndex, momDir)*grid.getDx(momDir)/2.0;
            integralResult[1] += coord3Val*baseIntegral;
            integralResult[2] += coord3Val*coord3Val*baseIntegral;
          }
          mom0MatrixVector[h](i, j) = integralResult[0];
          mom1MatrixVector[h](i, j) = integralResult[1];
          mom2MatrixVector[h](i, j) = integralResult[2];
        }
      }
    }

    // Get 3D Mass Matrix
    Lucee::Matrix<double> tempMassMatrix3d(nlocal3d, nlocal3d);
    nodalBasis3d->getMassMatrix(tempMassMatrix3d);
    Eigen::MatrixXd massMatrix3d(nlocal3d, nlocal3d);
    copyLuceeToEigen(tempMassMatrix3d, massMatrix3d);

    // Multiply matrices by inverse of mass matrix
    for (int h = 0; h < nlocal3d; h++)
    {
      mom0MatrixVector[h] = massMatrix3d.inverse()*mom0MatrixVector[h];
      mom1MatrixVector[h] = massMatrix3d.inverse()*mom1MatrixVector[h];
      mom2MatrixVector[h] = massMatrix3d.inverse()*mom2MatrixVector[h];
    }
  }

  Lucee::UpdaterStatus
  DistFuncMomentCalcWeighted3D::update(double t)
  {
    
    // get hold of grid
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // get input field (5d)
    const Lucee::Field<5, double>& distF = this->getInp<Lucee::Field<5, double> >(0);
    // get weighting field (3d)
    const Lucee::Field<3, double>& weightF = this->getInp<Lucee::Field<3, double> >(1);
    // get output field (3d)
    Lucee::Field<3, double>& momentOut = this->getOut<Lucee::Field<3, double> >(0);

    // create duplicate to store local moments
    Lucee::Field<3, double> moment = momentOut.duplicate();

    // local region to update (This is the 5D region. The 3D region is
    // assumed to have the same cell layout as the X-direction of the 5D region)
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();
    Lucee::Region<5, int> localExtRgn = distF.getExtRegion();

    // Make sure we integrate over conf. space ghost cells
    localRgn.setLower(0, localExtRgn.getLower(0));
    localRgn.setLower(1, localExtRgn.getLower(1));
    localRgn.setLower(2, localExtRgn.getLower(2));
    localRgn.setUpper(0, localExtRgn.getUpper(0));
    localRgn.setUpper(1, localExtRgn.getUpper(1));
    localRgn.setUpper(2, localExtRgn.getUpper(2));

    // clear out contents of output field
    moment = 0.0;

    // iterators into fields
    Lucee::ConstFieldPtr<double> distFPtr = distF.createConstPtr();
    Lucee::ConstFieldPtr<double> weightFPtr = weightF.createConstPtr();
    Lucee::FieldPtr<double> momentPtr = moment.createPtr();

    int idx[5];
    double xc[5];
    Lucee::RowMajorSequencer<5> seq(localRgn);
    unsigned nlocal3d = nodalBasis3d->getNumNodes();
    unsigned nlocal5d = nodalBasis5d->getNumNodes();

    int lower3d[] = {localRgn.getLower(0), localRgn.getLower(1), localRgn.getLower(2)};
    int upper3d[] = {localRgn.getUpper(0), localRgn.getUpper(1), localRgn.getUpper(2)};
    Lucee::Region<3, int> rgn3d(lower3d, upper3d);
    Lucee::RowMajorIndexer<3> idxr(rgn3d);
    Lucee::RowMajorSequencer<3> seq3d(rgn3d);

    int localPositionCells = localRgn.getShape(0)*localRgn.getShape(1)*localRgn.getShape(2);
    std::vector<double> localMoment(localPositionCells*nlocal3d);

    while(seq.step())
    {
      seq.fillWithIndex(idx);
      grid.setIndex(idx);
      grid.getCentroid(xc);

      moment.setPtr(momentPtr, idx[0], idx[1], idx[2]);
      distF.setPtr(distFPtr, idx);
      weightF.setPtr(weightFPtr, idx[0], idx[1], idx[2]);

      Eigen::VectorXd distfVec(nlocal5d);
      for (int i = 0; i < nlocal5d; i++)
        distfVec(i) = distFPtr[i];

      // Loop over each component of the weighting function
      for (int h = 0; h < nlocal3d; h++)
      {
        Eigen::VectorXd resultVector(nlocal3d);
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
        for (int i = 0; i < nlocal3d; i++)
          momentPtr[i] = momentPtr[i] + resultVector(i);
      }
    }

    int idx3d[3];
    // Loop over each 'local' position space cell
    while(seq3d.step())
    {
      seq3d.fillWithIndex(idx3d);
      int cellIndex = idxr.getIndex(idx3d);

      moment.setPtr(momentPtr, idx3d);
      // copy data to vector
      for (int i = 0; i < nlocal3d; i++)
        localMoment[cellIndex*nlocal3d+i] = momentPtr[i];
    }

    // Above loop computes moments on local phase-space domain. We need to
    // sum across velocity space to get total moment on configuration
    // space.
    std::vector<double> reducedMoment(localPositionCells*nlocal3d);
    // we need to get moment communicator of field as updater's moment
    // communicator is same as its grid's moment communicator. In this
    // case, grid is phase-space grid, which is not what we want.
    TxCommBase *momComm = momentOut.getMomComm();
    // amount to communicate
    momComm->allreduce(localPositionCells*nlocal3d, localMoment, reducedMoment, TX_SUM);

    seq3d.reset();
    Lucee::FieldPtr<double> momentOutPtr = momentOut.createPtr();
    // Copy reducedMoment to output field
    while(seq3d.step())
    {
      seq3d.fillWithIndex(idx3d);
      int cellIndex = idxr.getIndex(idx3d);
      momentOut.setPtr(momentOutPtr, idx3d);
      // copy data to vector
      for (int i = 0; i < nlocal3d; i++)
        momentOutPtr[i] = reducedMoment[cellIndex*nlocal3d+i];
    }

    return Lucee::UpdaterStatus();
  }

  void
  DistFuncMomentCalcWeighted3D::declareTypes()
  {
    // distribution function (5d)
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // weighting function (3d)
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // output field (3d)
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
  }

  void
  DistFuncMomentCalcWeighted3D::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}
