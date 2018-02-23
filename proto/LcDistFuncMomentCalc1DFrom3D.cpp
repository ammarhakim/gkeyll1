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
    momentLocal = 0;
  }

  DistFuncMomentCalc1DFrom3D::~DistFuncMomentCalc1DFrom3D()
  {
    delete momentLocal;
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
    Eigen::MatrixXd volCoords1d(nVolQuad1d, (unsigned) NC1);
    copyLuceeToEigen(tempVolCoords1d, volCoords1d);

    // get volume interpolation matrices for 3d element
    int nVolQuad3d = nodalBasis3d->getNumGaussNodes();
    std::vector<double> volWeights3d(nVolQuad3d);
    Lucee::Matrix<double> tempVolQuad3d(nVolQuad3d, nlocal3d);
    Lucee::Matrix<double> tempVolCoords3d(nVolQuad3d, (unsigned) NC3);

    nodalBasis3d->getGaussQuadData(tempVolQuad3d, tempVolCoords3d, volWeights3d);

    Eigen::MatrixXd volQuad3d(nVolQuad3d, nlocal3d);
    copyLuceeToEigen(tempVolQuad3d, volQuad3d);
    // TESTING STUFF
    Eigen::MatrixXd volCoords3d(nVolQuad3d, (unsigned) NC3);
    copyLuceeToEigen(tempVolCoords3d, volCoords3d);

    mom0Matrix = Eigen::MatrixXd::Zero(nlocal1d, nlocal3d);
    mom1Matrix = Eigen::MatrixXd::Zero(nlocal1d, nlocal3d);
    mom2Matrix = Eigen::MatrixXd::Zero(nlocal1d, nlocal3d);

    for (int i = 0; i < nlocal1d; i++)
    {
      for (int j = 0; j < nlocal3d; j++)
      {
        // Compute integral of phi1d_i * phi3d_j
        double integralResult[3] = {};
        for (int gaussIndex = 0; gaussIndex < volWeights3d.size(); gaussIndex++)
        {
          // volWeights is created with a factor of grid.getDx(dim). Here we divide
          // by this factor along the non-uniform dimension and account for it when
          // looping over cells in the update part below. 
          double baseIntegral = volWeights3d[gaussIndex]*volQuad1d(gaussIndex % nVolQuad1d, i)*
            volQuad3d(gaussIndex, j)/grid.getDx(2);
          integralResult[0] += baseIntegral;
          // Get coordinate of quadrature point in direction momDir.
          // Will multiply by Dx in update section to support non-uniform grid.
          double coord2Val = volCoords3d(gaussIndex, momDir)*0.50; //grid.getDx(momDir);
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
    Lucee::Field<1, double>& momentGlobal = this->getOut<Lucee::Field<1, double> >(0);

    // clear out contents of output field
    if (!momentLocal)
    {
      // allocate memory for local moment calculation if not already
      // done: we need to ensure space is also allocated for the
      // ghost-cells as otherwise there is a size mis-match in the
      // allReduce call to sync across velocity space
      Lucee::Region<1, int> localRgn = momentGlobal.getRegion();
      Lucee::Region<1, int> localExtRgn = momentGlobal.getExtRegion();
      
      int lowerConf[1];
      int upperConf[1];
      int lg[1];
      int ug[1];
      for (int i=0; i < 1; ++i)
      {
        lowerConf[i] = localRgn.getLower(i);
        upperConf[i] = localRgn.getUpper(i);
        lg[i] = localRgn.getLower(i) - localExtRgn.getLower(i);
        ug[i] = localExtRgn.getUpper(i) - localRgn.getUpper(i);
      }
      Lucee::Region<1, int> rgnConf(lowerConf, upperConf);
      momentLocal = new Lucee::Field<1, double>(rgnConf, momentGlobal.getNumComponents(), lg, ug);
    }

    // local region to update (This is the 3D region. The 1D region is
    // assumed to have the same cell layout as the X-direction of the 3D region)
    Lucee::Region<3, int> localRgn = grid.getLocalRegion();

    // clear out contents of output field
    (*momentLocal) = 0.0;

    // iterators into fields
    Lucee::ConstFieldPtr<double> distFPtr = distF.createConstPtr();
    Lucee::FieldPtr<double> momentPtr = momentLocal->createPtr();

    int idx[3];
    double xc[3];
    Lucee::RowMajorSequencer<3> seq(localRgn);
    unsigned nlocal1d = nodalBasis1d->getNumNodes();
    unsigned nlocal3d = nodalBasis3d->getNumNodes();

    int localPositionCells = momentGlobal.getExtRegion().getVolume();

    while(seq.step())
    {
      seq.fillWithIndex(idx);

      grid.setIndex(idx);
      grid.getCentroid(xc);
      
      double physicalVol = grid.getVolume();	// Volume of this cell (physical units).
      double DxMu        = physicalVol/grid.getSurfArea(2); // Mu cell length.
      double momDirDx    = physicalVol/grid.getSurfArea(momDir); // Cell length along momDir.

      momentLocal->setPtr(momentPtr, idx[0]);
      distF.setPtr(distFPtr, idx);

      Eigen::VectorXd distfVec(nlocal3d);
      for (int i = 0; i < nlocal3d; i++)
        distfVec(i) = distFPtr[i];

      // Accumulate contribution to momentGlobal from this cell
      Eigen::VectorXd resultVector(nlocal1d);

      // Here volume/SurfArea gives the correct cell length along the non-uniform dimension.
      if (calcMom == 0)
        resultVector.noalias() = mom0Matrix*distfVec*DxMu;
      else if (calcMom == 1)
        resultVector.noalias() = (mom1Matrix*momDirDx + xc[momDir]*mom0Matrix)*DxMu*distfVec;
      else if (calcMom == 2)
        resultVector.noalias() = (mom2Matrix*momDirDx*momDirDx + 2.0*xc[momDir]*mom1Matrix*momDirDx +
          xc[momDir]*xc[momDir]*mom0Matrix)*DxMu*distfVec;

      for (int i = 0; i < nlocal1d; i++)
        momentPtr[i] = momentPtr[i] + resultVector(i);
    }

    // Above loop computes moments on local phase-space domain. We need to
    // sum across velocity space to get total momentLocal on configuration
    // space.
    
    // we need to get moment communicator of field as updater's moment
    // communicator is same as its grid's moment communicator. In this
    // case, grid is phase-space grid, which is not what we want.
    TxCommBase *momComm = momentGlobal.getMomComm();
    unsigned xsize = localPositionCells*nlocal1d; // amount to communicate
    momComm->allreduce(xsize, &momentLocal->first(), &momentGlobal.first(), TX_SUM);

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
