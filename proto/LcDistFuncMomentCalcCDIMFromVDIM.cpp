/**
 * @file	LcDistFuncMomentCalcCDIMFromVDIM.cpp
 *
 * @brief	Updater to compute arbitrary configuration space moments of an arbitrary phase space distribution function
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDistFuncMomentCalcCDIMFromVDIM.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <cmath>
#include <vector>

namespace Lucee
{
// set id for module system
  template <> const char *DistFuncMomentCalcCDIMFromVDIM<1,1>::id = "DistFuncMomentCalc1X1V";
  template <> const char *DistFuncMomentCalcCDIMFromVDIM<1,2>::id = "DistFuncMomentCalc1X2V";
  template <> const char *DistFuncMomentCalcCDIMFromVDIM<1,3>::id = "DistFuncMomentCalc1X3V";
  template <> const char *DistFuncMomentCalcCDIMFromVDIM<2,2>::id = "DistFuncMomentCalc2X2V";
  template <> const char *DistFuncMomentCalcCDIMFromVDIM<2,3>::id = "DistFuncMomentCalc2X3V";
  //template <> const char *DistFuncMomentCalcCDIMFromVDIM<3,3>::id = "DistFuncMomentCalc3X3V";

// makes indexing a little more sane
  static const unsigned IX = 0;
  static const unsigned IY = 1;
  static const unsigned IZ = 2;

  template <unsigned CDIM, unsigned VDIM>
  DistFuncMomentCalcCDIMFromVDIM<CDIM, VDIM>::DistFuncMomentCalcCDIMFromVDIM()
    : Lucee::UpdaterIfc()
  {
  }

  template <unsigned CDIM, unsigned VDIM>
  DistFuncMomentCalcCDIMFromVDIM<CDIM, VDIM>::~DistFuncMomentCalcCDIMFromVDIM()
  {
    delete momentLocal;
  }

  template <unsigned CDIM, unsigned VDIM>  
  void
  DistFuncMomentCalcCDIMFromVDIM<CDIM, VDIM>::readInput(Lucee::LuaTable& tbl)
  {
    const unsigned NDIM = CDIM+VDIM;
    // call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    // get hold of phase space element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("phaseBasis"))
      phaseBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("phaseBasis");
    else
      throw Lucee::Except(
        "DistFuncMomentCalcCDIMFromVDIM::readInput: Must specify phase-space basis using 'phaseBasis'");

    // get hold of configuration space element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<CDIM> >("confBasis"))
      confBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<CDIM> >("confBasis");
    else
      throw Lucee::Except(
        "DistFuncMomentCalcCDIMFromVDIM::readInput: Must specify configuration-space basis using 'confBasis'");

    // get moment to compute
    if (tbl.hasNumber("moment"))
    calcMom = (unsigned) tbl.getNumber("moment");
    else
      throw Lucee::Except(
        "DistFuncMomentCalcCDIMFromVDIM::readInput: Must specify moment using 'moment'");

    if (calcMom > 2)
    {
      Lucee::Except lce("DistFuncMomentCalcCDIMFromVDIM::readInput: Only 'moment' 0, 1, or 2 is supported. ");
      lce << "Supplied " << calcMom << " instead";
      throw lce;
    }
  }

  template <unsigned CDIM, unsigned VDIM>
  void
  DistFuncMomentCalcCDIMFromVDIM<CDIM, VDIM>::initialize()
  {
    const unsigned NDIM = CDIM+VDIM;
    // call base class method
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    // local region to update
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    // get number of nodes in configuration space and phase space
    unsigned nlocalConf = confBasis->getNumNodes();
    unsigned nlocalPhase = phaseBasis->getNumNodes();

    // get volume interpolation matrices for configuration space element
    int nVolQuadConf = confBasis->getNumGaussNodes();
    std::vector<double> volWeightsConf(nVolQuadConf);
    Lucee::Matrix<double> tempVolQuadConf(nVolQuadConf, nlocalConf);
    Lucee::Matrix<double> tempVolCoordsConf(nVolQuadConf, CNC);

    confBasis->getGaussQuadData(tempVolQuadConf, tempVolCoordsConf, volWeightsConf);

    Eigen::MatrixXd volQuadConf(nVolQuadConf, nlocalConf);
    copyLuceeToEigen(tempVolQuadConf, volQuadConf);
    // TESTING STUFF
    Eigen::MatrixXd volCoordsConf(nVolQuadConf, (unsigned) CNC);
    copyLuceeToEigen(tempVolCoordsConf, volCoordsConf);

    // get volume interpolation matrices for phase space element
    int nVolQuadPhase = phaseBasis->getNumGaussNodes();
    std::vector<double> volWeightsPhase(nVolQuadPhase);
    Lucee::Matrix<double> tempVolQuadPhase(nVolQuadPhase, nlocalPhase);
    Lucee::Matrix<double> tempVolCoordsPhase(nVolQuadPhase, (unsigned) PNC);

    phaseBasis->getGaussQuadData(tempVolQuadPhase, tempVolCoordsPhase, volWeightsPhase);

    Eigen::MatrixXd volQuadPhase(nVolQuadPhase, nlocalPhase);
    copyLuceeToEigen(tempVolQuadPhase, volQuadPhase);
    // TESTING STUFF
    Eigen::MatrixXd volCoordsPhase(nVolQuadPhase, (unsigned) PNC);
    copyLuceeToEigen(tempVolCoordsPhase, volCoordsPhase);

    // Get configuration space Mass Matrix
    Lucee::Matrix<double> tempMassMatrixConf(nlocalConf, nlocalConf);
    confBasis->getMassMatrix(tempMassMatrixConf);
    Eigen::MatrixXd massMatrixConf(nlocalConf, nlocalConf);
    copyLuceeToEigen(tempMassMatrixConf, massMatrixConf);

    nMom = 1;
    if (calcMom == 1)
      nMom = VDIM;
    else if (calcMom == 2)
      nMom = VDIM*(VDIM+1)/2;

    mom0Matrix = Eigen::MatrixXd::Zero(nlocalConf, nlocalPhase);
    // matrices needed to compute the first moment
    mom1Matrix = std::vector<Eigen::MatrixXd>(VDIM);
    // matrices needed to compute the second moment
    mom2Matrix = std::vector< std::vector<Eigen::MatrixXd> >(VDIM*(VDIM+1)/2);
    // momDir is direction of moment, 1 = X, 2 = Y, 3 = Z
    int momDir;
    // momDir2 is direction of other moment, 1 = X, 2 = Y, 3 = Z, for computing second rank tensor components
    int momDir2;
    double integralResultZerothMoment;
    double integralResultFirstMoment[VDIM];
    // number of components in a symmetric second rank tensor
    double integralResultSecondMoment[VDIM*(VDIM+1)/2];
    // counter needed for indexing since tensor is symmetric
    int ctr = 0;

    for (int i = 0; i < nlocalConf; ++i)
    {
      for (int j = 0; j < nlocalPhase; ++j)
      {
        integralResultZerothMoment = 0;
        // Compute integral of phiConf_i * phiPhase_j
        for (int gaussIndex = 0; gaussIndex < volWeightsPhase.size(); ++gaussIndex)
        {
          double baseIntegral = volWeightsPhase[gaussIndex]*volQuadConf(gaussIndex % nVolQuadConf, i)*
          volQuadPhase(gaussIndex, j);
          integralResultZerothMoment += baseIntegral;
        }
        mom0Matrix(i, j) = integralResultZerothMoment;
      }
    }
    // Multiply matrices by inverse of mass matrix
    mom0Matrix = massMatrixConf.inverse()*mom0Matrix;

    for (int h = 0; h < nMom; ++h)
    {
      mom1Matrix[h] = Eigen::MatrixXd::Zero(nlocalConf, nlocalPhase);
      momDir = h+1;
      // Compute integral of phiConf_i * phiPhase_j
      for (int i = 0; i < nlocalConf; ++i)
      {
        for (int j = 0; j < nlocalPhase; ++j)
        {
          integralResultFirstMoment[h] = 0.0;
          for (int gaussIndex = 0; gaussIndex < volWeightsPhase.size(); ++gaussIndex)
          {
            double baseIntegral = volWeightsPhase[gaussIndex]*volQuadConf(gaussIndex % nVolQuadConf, i)*
              volQuadPhase(gaussIndex, j);
            // Get coordinate of quadrature point in direction momDir
            double coord2Val = volCoordsPhase(gaussIndex, momDir)*grid.getDx(momDir)/2.0;
            integralResultFirstMoment[h] += coord2Val*baseIntegral;
          }
          mom1Matrix[h](i, j) = integralResultFirstMoment[h];
        }
      }
      // Multiply matrices by inverse of mass matrix
      mom1Matrix[h] = massMatrixConf.inverse()*mom1Matrix[h];
    }

    for (int h = 0; h < VDIM; h++)
    {
      for (int g = h; g < VDIM; g++)
      {
      mom2Matrix[h][g] = Eigen::MatrixXd::Zero(nlocalConf, nlocalPhase);
      momDir = g+1;
      momDir2 = h+1;
      // Compute integral of phiConf_i * phiPhase_j
      for (int i = 0; i < nlocalConf; i++)
      {
        for (int j = 0; j < nlocalPhase; j++)
        {
          for (int gaussIndex = 0; gaussIndex < volWeightsPhase.size(); gaussIndex++)
          {
            integralResultSecondMoment[ctr] = 0.0;
            double baseIntegral = volWeightsPhase[gaussIndex]*volQuadConf(gaussIndex % nVolQuadConf, i)*
              volQuadPhase(gaussIndex, j);
            // Get coordinate of quadrature point in direction momDir
            double coord2Val1 = volCoordsPhase(gaussIndex, momDir)*grid.getDx(momDir)/2.0;
            double coord2Val2 = volCoordsPhase(gaussIndex, momDir2)*grid.getDx(momDir2)/2.0;
            integralResultSecondMoment[ctr] += coord2Val1*coord2Val2*baseIntegral;
          }
          mom2Matrix[h][g](i, j) = integralResultSecondMoment[ctr];
        }
      }
      ctr = ctr + 1;
      // Multiply matrices by inverse of mass matrix
      mom2Matrix[h][g] = massMatrixConf.inverse()*mom2Matrix[h][g];
      }
    }

    int lowerConf[CDIM];
    int upperConf[CDIM];
    int lg[CDIM];
    int ug[CDIM];
    for (int i=0; i<CDIM; ++i)
    {
      lowerConf[i] = localRgn.getLower(i);
      upperConf[i] = localRgn.getUpper(i);
      lg[i] = 0;
      ug[i] = 0;
    }

    Lucee::Region<CDIM, int> rgnConf(lowerConf, upperConf);

    // allocate space for storing local moment calculation
    momentLocal = new Lucee::Field<CDIM, double>(rgnConf, nlocalConf*nMom, lg, ug);
  }

  template <unsigned CDIM, unsigned VDIM>
  Lucee::UpdaterStatus
  DistFuncMomentCalcCDIMFromVDIM<CDIM, VDIM>::update(double t)
  {
    const unsigned NDIM = CDIM+VDIM;
    // get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    // get input field (CDIM + VDIM)
    const Lucee::Field<NDIM, double>& distF = this->getInp<Lucee::Field<NDIM, double> >(0);
    // get output field (CDIM)
    Lucee::Field<CDIM, double>& momentGlobal = this->getOut<Lucee::Field<CDIM, double> >(0);

    // local region to update (This is the CDIM + VDIM region. The CDIM region is
    // assumed to have the same cell layout as the X-direction of the CDIM+VDIM region)
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    int idx[NDIM];
    double xc[NDIM];
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    unsigned nlocalConf = confBasis->getNumNodes();
    unsigned nlocalPhase = phaseBasis->getNumNodes();

    // calculate total number of local configuration space cells
    int localPositionCells = 1;
    for (unsigned d=0; d<CDIM; ++d)
      localPositionCells *= localRgn.getShape(d);

    // clear out contents of output field
    (*momentLocal) = 0.0;

    // iterators into fields
    Lucee::ConstFieldPtr<double> distFPtr = distF.createConstPtr();
    Lucee::FieldPtr<double> momentPtr = momentLocal->createPtr();
    Eigen::VectorXd distfVec(nlocalPhase);

    // Accumulate contribution to moment from this cell
    
    std::vector<Eigen::VectorXd> resultVector = std::vector<Eigen::VectorXd>(nMom);
    for (int i = 0; i< nMom; ++i)
    {
      resultVector[i] = Eigen::VectorXd::Zero(nlocalConf);
    }
    while(seq.step())
    {
      seq.fillWithIndex(idx);

      grid.setIndex(idx);
      grid.getCentroid(xc);

      momentLocal->setPtr(momentPtr, idx);
      distF.setPtr(distFPtr, idx);

      for (int i = 0; i < nlocalPhase; i++)
        distfVec(i) = distFPtr[i];


      if (calcMom == 0)
      {
        resultVector[0].noalias() = mom0Matrix*distfVec;
      }
      else if (calcMom == 1)
      {
        int momDir;
        for (int h = 0; h < nMom; h++)
        {
          momDir = h+1;
          resultVector[h].noalias() = (mom1Matrix[h] + xc[momDir]*mom0Matrix)*distfVec;
        }
      }
      else if (calcMom == 2)
      {
        int momDir;
        int momDir2;
        int ctr = 0;
        for (int h = 0; h < VDIM; h++)
        {
          for (int g = h; g < VDIM; g++)
          {
            momDir = g+1;
            momDir2 = h+1;
            resultVector[ctr++].noalias() = (mom2Matrix[h][g] + 2*xc[momDir]*mom1Matrix[h] + xc[momDir]*xc[momDir2]*mom0Matrix)*distfVec;
          }
        }
      }

      for (int i = 0; i < nlocalConf; i++)
        for (int j = 0; j < nMom; ++j)
          momentPtr[i*nMom+j] += resultVector[j](i);
    }

// Above loop computes moments on local phase-space domain. We need to
// sum across velocity space to get total moment on configuration
// space.

// we need to get moment communicator of field as updater's moment
// communicator is same as its grid's moment communicator. In this
// case, grid is phase-space grid, which is not what we want.
    TxCommBase *momComm = momentGlobal.getMomComm();
    unsigned xsize = localPositionCells*nlocalConf*nMom; // amount to communicate
    momComm->allreduce(xsize, &momentLocal->firstInterior(), &momentGlobal.firstInterior(), TX_SUM);

    return Lucee::UpdaterStatus();
  }

  template <unsigned CDIM, unsigned VDIM>  
  void
  DistFuncMomentCalcCDIMFromVDIM<CDIM, VDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<CDIM+VDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<CDIM, double>));
  }

  template <unsigned CDIM, unsigned VDIM>
  void
  DistFuncMomentCalcCDIMFromVDIM<CDIM, VDIM>::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }

// instantiations
  template class DistFuncMomentCalcCDIMFromVDIM<1,1>;
  template class DistFuncMomentCalcCDIMFromVDIM<1,2>;
  template class DistFuncMomentCalcCDIMFromVDIM<1,3>;
  template class DistFuncMomentCalcCDIMFromVDIM<2,2>;
  template class DistFuncMomentCalcCDIMFromVDIM<2,3>;  
}
