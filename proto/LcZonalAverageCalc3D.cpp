/**
 * @file	LcZonalAverageCalc3D.cpp
 *
 * @brief	Updater to compute zonal (y) average of a 3d potential.
 * Output will be a 2d field.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcZonalAverageCalc3D.h>
#include <LcGlobals.h>
#include <LcLinAlgebra.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <cmath>
#include <vector>

namespace Lucee
{
  const char *ZonalAverageCalc3D::id = "ZonalAverageCalc3D";

  ZonalAverageCalc3D::ZonalAverageCalc3D()
    : Lucee::UpdaterIfc()
  {
  }

  ZonalAverageCalc3D::~ZonalAverageCalc3D()
  {
    delete moment;
  }
  
  void
  ZonalAverageCalc3D::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    // get hold of 3d element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except(
        "ZonalAverageCalc3D::readInput: Must specify 3D element to use using 'basis3d'");

    // get hold of 2D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis2d"))
      nodalBasis2d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis2d");
    else
      throw Lucee::Except("ZonalAverageCalc3D::readInput: Must specify 2D element to use using 'basis2d'");
  }

  void
  ZonalAverageCalc3D::initialize()
  {
    // call base class method
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    // local region to update
    Lucee::Region<3, int> localRgn = grid.getLocalRegion();

    // get number of nodes in 1D and 2D
    unsigned nlocal2d = nodalBasis2d->getNumNodes();
    unsigned nlocal3d = nodalBasis3d->getNumNodes();

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

    // get volume interpolation matrices for 3d element
    int nVolQuad3d = nodalBasis3d->getNumGaussNodes();
    std::vector<double> volWeights3d(nVolQuad3d);
    Lucee::Matrix<double> tempVolQuad3d(nVolQuad3d, nlocal3d);
    Lucee::Matrix<double> tempVolCoords3d(nVolQuad3d, 3);

    nodalBasis3d->getGaussQuadData(tempVolQuad3d, tempVolCoords3d, volWeights3d);

    Eigen::MatrixXd volQuad3d(nVolQuad3d, nlocal3d);
    copyLuceeToEigen(tempVolQuad3d, volQuad3d);
    Eigen::MatrixXd volCoords3d(nVolQuad3d, 3);
    copyLuceeToEigen(tempVolCoords3d, volCoords3d);

    mom0Matrix = Eigen::MatrixXd(nlocal2d, nlocal3d);

    for (int i = 0; i < nlocal2d; i++)
    {
      for (int j = 0; j < nlocal3d; j++)
      {
        double integralResult = 0.0;
        // Compute integral of phi2d_i * phi3d_j in 3d
        for (int gaussIndex3d = 0; gaussIndex3d < volWeights3d.size(); gaussIndex3d++)
        {
          bool quadMatchFound = false;
          int quadMatchIndex = 0;
          // Need to find correct match to this gaussIndex in 2d
          for (int gaussIndex2d = 0; gaussIndex2d < volWeights2d.size(); gaussIndex2d++)
          {
            if (std::fabs(volCoords2d(gaussIndex2d,0)-volCoords3d(gaussIndex3d,0)) < 1e-10 &&
                std::fabs(volCoords2d(gaussIndex2d,1)-volCoords3d(gaussIndex3d,2)) < 1e-10)
            {
              quadMatchFound = true;
              quadMatchIndex = gaussIndex2d;
            }
          }
          if (quadMatchFound == true)
            integralResult += volWeights3d[gaussIndex3d]*volQuad2d(quadMatchIndex, i)*volQuad3d(gaussIndex3d, j);
          else
            throw Lucee::Except("ZonalAverageCalc3D::initialize: Unable to find a corresponding 2d quadrature point.");
        }
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

    int lower[2], upper[2];
    lower[0] = localRgn.getLower(0); upper[0] = localRgn.getUpper(0);
    lower[1] = localRgn.getLower(2); upper[1] = localRgn.getUpper(2);
    Lucee::Region<2, int> local2D(lower, upper);
    int lg[2] = {0,0}, ug[2] = {0,0};
    // allocate space for storing local moment calculation
    moment = new Lucee::Field<2, double>(local2D, nlocal2d, lg, ug);
  }

  Lucee::UpdaterStatus
  ZonalAverageCalc3D::update(double t)
  {
    
    // get hold of grid
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    // get input field (3d)
    const Lucee::Field<3, double>& phiIn = this->getInp<Lucee::Field<3, double> >(0);
    // get output field (2d)
    Lucee::Field<2, double>& phiZonalOut = this->getOut<Lucee::Field<2, double> >(0);

  // local region to update (This is the 3d region. The 2D region is
  // assumed to have the same cell layout as the X-direction of the 3d region)
    Lucee::Region<3, int> localRgn = grid.getLocalRegion();

    // clear out contents of output field
    (*moment) = 0.0;

    // iterators into fields
    Lucee::ConstFieldPtr<double> phiPtr = phiIn.createConstPtr();
    Lucee::FieldPtr<double> momentPtr = moment->createPtr();
 
    int idx[3];
    Lucee::RowMajorSequencer<3> seq(localRgn);
    unsigned nlocal2d = nodalBasis2d->getNumNodes();
    unsigned nlocal3d = nodalBasis3d->getNumNodes();

    while(seq.step())
    {
      seq.fillWithIndex(idx);

      grid.setIndex(idx);
      // phi at (ix,iy,iz) gets accumulated to phiAvg(ix,iz)
      moment->setPtr(momentPtr, idx[0], idx[2]);
      phiIn.setPtr(phiPtr, idx);

      Eigen::VectorXd phiVec(nlocal3d);
      for (int i = 0; i < nlocal3d; i++)
        phiVec(i) = phiPtr[i];

      // Accumulate contribution to moment from this cell
      Eigen::VectorXd resultVector = mom0Matrix*phiVec;

      for (int i = 0; i < nlocal2d; i++)
        momentPtr[i] = momentPtr[i] + resultVector(i);
    }

    int lower2d[] = {localRgn.getLower(0), localRgn.getLower(2)};
    int upper2d[] = {localRgn.getUpper(0), localRgn.getUpper(2)};
    Lucee::Region<2, int> rgn2d(lower2d, upper2d);
    Lucee::RowMajorIndexer<2> idxr(rgn2d);
    Lucee::RowMajorSequencer<2> seq2d(rgn2d);
    int idx2d[2];
    int localPositionCells = localRgn.getShape(0)*localRgn.getShape(2);
    std::vector<double> localMoment(localPositionCells*nlocal2d);
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
 
    // we need to get moment communicator of field as updater's moment
    // communicator is same as its grid's moment communicator. In this
    // case, grid is phase-space grid, which is not what we want.
    TxCommBase *momComm = phiZonalOut.getMomComm();
    // amount to communicate
    unsigned datasize = phiZonalOut.getNumComponents()*localRgn.getShape(0)*localRgn.getShape(2);
    std::vector<double> reducedMoment(localPositionCells*nlocal2d);
    momComm->allreduce(datasize, localMoment, reducedMoment, TX_SUM);
    //momComm->allreduce(datasize, &moment->firstInterior(), &phiZonalOut.firstInterior(), TX_SUM);

    seq2d.reset();
    Lucee::FieldPtr<double> phiZonalOutPtr = phiZonalOut.createPtr();
    // Copy reducedMoment to output field
    while(seq2d.step())
    {
      seq2d.fillWithIndex(idx2d);
      int cellIndex = idxr.getIndex(idx2d);
      phiZonalOut.setPtr(phiZonalOutPtr, idx2d);
      // copy data to vector
      for (int i = 0; i < nlocal2d; i++)
        phiZonalOutPtr[i] = reducedMoment[cellIndex*nlocal2d+i];
    }

    return Lucee::UpdaterStatus();
  }

  void
  ZonalAverageCalc3D::declareTypes()
  {
    // Input potential on a 3d (x,y,z) field
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Zonal-averaged potential on a 2d (x,z) field
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }

  void
  ZonalAverageCalc3D::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}
