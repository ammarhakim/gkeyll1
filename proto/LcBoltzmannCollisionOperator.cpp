/**
 * @file	LcBoltzmannCollisionOperator.cpp
 *
 * @brief Updater to compute a collision operator for RHS of the Bolztmann 
 *  as defined in Meier & Shumlak (2012) [http://dx.doi.org/10.1063/1.4736975]
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBoltzmannCollisionOperator.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <cmath>
#include <vector>

namespace Lucee
{
// set id for module system
  template <> const char *BoltzmannCollisionOperator<1,1>::id =
    "BoltzmannCollisionOperator1X1V";
  template <> const char *BoltzmannCollisionOperator<1,2>::id = 
    "BoltzmannCollisionOperator1X2V";
  template <> const char *BoltzmannCollisionOperator<1,3>::id = 
    "BoltzmannCollisionOperator1X3V";
  template <> const char *BoltzmannCollisionOperator<2,2>::id = 
    "BoltzmannCollisionOperator2X2V";
  template <> const char *BoltzmannCollisionOperator<2,3>::id = 
    "BoltzmannCollisionOperator2X3V";
  //template <> const char *DistFuncMomentCalcCDIMFromVDIM<3,3>::id = "DistFuncMomentCalc3X3V";


  template <unsigned CDIM, unsigned VDIM>
  bool
  BoltzmannCollisionOperator<CDIM,VDIM>::sameConfigCoords(unsigned n, unsigned cn, double dxMin,
    const Eigen::MatrixXd& phaseC, const Eigen::MatrixXd& confC)
  {
    for (unsigned dim = 0; dim<CDIM; ++dim)
      if (! (std::fabs(phaseC(n, dim)-confC(cn, dim))<1e-4*dxMin) )
        return false;
    return true;
  }

  template <unsigned CDIM, unsigned VDIM>
  BoltzmannCollisionOperator<CDIM, VDIM>::BoltzmannCollisionOperator()
    : Lucee::UpdaterIfc()
  {
    momentLocal = 0;
  }

  template <unsigned CDIM, unsigned VDIM>
  BoltzmannCollisionOperator<CDIM, VDIM>::~BoltzmannCollisionOperator()
  {
    delete momentLocal;
  }

  template <unsigned CDIM, unsigned VDIM>  
  void
  BoltzmannCollisionOperator<CDIM, VDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);    
    const unsigned NDIM = CDIM+VDIM;

    // get hold of phase space element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("phaseBasis"))
      phaseBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("phaseBasis");
    else
      throw Lucee::Except(
        "BoltzmannCollisionOperator::readInput: Must specify phase-space basis using 'phaseBasis'");

    // get hold of configuration space element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<CDIM> >("confBasis"))
      confBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<CDIM> >("confBasis");
    else
      throw Lucee::Except(
        "BoltzmannCollisionOperator::readInput: Must specify configuration-space basis using 'confBasis'");

    // get moment to compute
    if (tbl.hasNumber("moment"))
    calcMom = (unsigned) tbl.getNumber("moment");
    else
      throw Lucee::Except(
        "BoltzmannCollisionOperator::readInput: Must specify moment using 'moment'");
  }

  template <unsigned CDIM, unsigned VDIM>
  void
  BoltzmannCollisionOperator<CDIM, VDIM>::initialize()
  {
    Lucee::UpdaterIfc::initialize();    
    const unsigned NDIM = CDIM+VDIM;

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
    Eigen::MatrixXd volCoordsPhase(nVolQuadPhase, (unsigned) PNC);
    copyLuceeToEigen(tempVolCoordsPhase, volCoordsPhase);

    // compute mapping of phase-space quadrature nodes to configuration space
    // quadrature nodes. The assumption here is that the node layout in phase-space
    // and configuration space are such that each node in phase-space has
    // exactly one node co-located with it in configuration space. No
    // "orphan" phase-space node are allowed, and an exception is thrown
    // if that occurs.
    phaseConfMap.resize(volWeightsPhase.size());

    double dxMin = grid.getDx(0);
    for (unsigned d=1; d<CDIM; ++d)
      dxMin = std::min(dxMin, grid.getDx(d));

    for (unsigned n=0; n<volWeightsPhase.size(); ++n)
    {
      bool pcFound = false;
      for (unsigned cn=0; cn<volWeightsConf.size(); ++cn)
        if (sameConfigCoords(n, cn, dxMin, volCoordsPhase, volCoordsConf))
        {
          phaseConfMap[n] = cn;
          pcFound = true;
          break;
        }
      if (!pcFound)
      {
        Lucee::Except lce(
          "BoltzmannCollisionOperator::initialize: No matching configuration space quadrature node for phase-space quadrature node ");
        lce << n;
        throw lce;
      }
    }

    // Get configuration space mass matrix
    Lucee::Matrix<double> tempMassMatrixConf(nlocalConf, nlocalConf);
    confBasis->getMassMatrix(tempMassMatrixConf);
    Eigen::MatrixXd massMatrixConf(nlocalConf, nlocalConf);
    copyLuceeToEigen(tempMassMatrixConf, massMatrixConf);

    nMom = 1;
    if (calcMom == 1)
      nMom = VDIM;
    nMom = 1;

    mom0Matrix = Eigen::MatrixXd::Zero(nlocalConf, nlocalPhase);    
    for (int i=0; i<nlocalConf; ++i)
    {
      for (int j=0; j<nlocalPhase; ++j)
      {
        double integralResultZerothMoment = 0;
        // Compute integral of phiConf_i * phiPhase_j
        for (int gaussIndex = 0; gaussIndex < volWeightsPhase.size(); ++gaussIndex)
          integralResultZerothMoment
            += volWeightsPhase[gaussIndex]*volQuadConf(phaseConfMap[gaussIndex],i)*volQuadPhase(gaussIndex,j);
        mom0Matrix(i, j) = integralResultZerothMoment;
      }
    }
    // Multiply matrices by inverse of mass matrix
    mom0Matrix = massMatrixConf.inverse()*mom0Matrix;

    mom1Matrix = std::vector<Eigen::MatrixXd>(VDIM);    
    for (int h=0; h<VDIM; ++h)
    {
      mom1Matrix[h] = Eigen::MatrixXd::Zero(nlocalConf, nlocalPhase);
      int momDir = h+CDIM;
      for (int i=0; i<nlocalConf; ++i)
      {
        for (int j=0; j<nlocalPhase; ++j)
        {
          double integralResultFirstMoment = 0.0;
          for (int gaussIndex = 0; gaussIndex < volWeightsPhase.size(); ++gaussIndex)
          {
            double baseIntegral
              = volWeightsPhase[gaussIndex]*volQuadConf(phaseConfMap[gaussIndex],i)*volQuadPhase(gaussIndex,j);
            // get coordinate of quadrature point in direction momDir
            double coord2Val = volCoordsPhase(gaussIndex, momDir)*grid.getDx(momDir)/2.0;
            integralResultFirstMoment += coord2Val*baseIntegral;
          }
          mom1Matrix[h](i,j) = integralResultFirstMoment;
        }
      }
      // multiply matrices by inverse of mass matrix
      mom1Matrix[h] = massMatrixConf.inverse()*mom1Matrix[h];
    }

    mom2Matrix = std::vector<Eigen::MatrixXd>(VDIM*(VDIM+1)/2);  
    int ctr = 0; // counter needed for indexing since tensor is symmetric    
    for (int h=0; h<VDIM; ++h)
    {
      for (int g=h; g<VDIM; ++g)
      {
        mom2Matrix[ctr] = Eigen::MatrixXd::Zero(nlocalConf, nlocalPhase);
        int momDir = h+CDIM;
        int momDir2 = g+CDIM;
        // Compute integral of phiConf_i * phiPhase_j
        for (int i=0; i<nlocalConf; ++i)
        {
          for (int j=0; j<nlocalPhase; ++j)
          {
            double integralResultSecondMoment = 0.0;
            for (int gaussIndex = 0; gaussIndex < volWeightsPhase.size(); ++gaussIndex)
            {
              double baseIntegral 
                = volWeightsPhase[gaussIndex]*volQuadConf(phaseConfMap[gaussIndex], i)*volQuadPhase(gaussIndex, j);
              // Get coordinate of quadrature point in direction momDir
              double coord2Val1 = volCoordsPhase(gaussIndex, momDir)*grid.getDx(momDir)/2.0;
              double coord2Val2 = volCoordsPhase(gaussIndex, momDir2)*grid.getDx(momDir2)/2.0;
              integralResultSecondMoment += coord2Val1*coord2Val2*baseIntegral;
            }
            mom2Matrix[ctr](i,j) = integralResultSecondMoment;
          }
        }
        // Multiply matrices by inverse of mass matrix
        mom2Matrix[ctr] = massMatrixConf.inverse()*mom2Matrix[ctr];
        ctr = ctr+1;
      }
    }
  }

  template <unsigned CDIM, unsigned VDIM>
  Lucee::UpdaterStatus
  BoltzmannCollisionOperator<CDIM, VDIM>::update(double t)
  {
    const unsigned NDIM = CDIM+VDIM;
    // get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    // get input field (CDIM + VDIM) 
    const Lucee::Field<NDIM, double>& distF = this->getInp<Lucee::Field<NDIM, double> >(0);
    // get output field (CDIM)
    Lucee::Field<CDIM, double>& momentGlobal = this->getOut<Lucee::Field<CDIM, double> >(0);

    if (!momentLocal)
    {
      // allocate memory for local moment calculation if not already
      // done: we need to ensure space is also allocated for the
      // ghost-cells as otherwise there is a size mis-match in the
      // allReduce call to sync across velocity space
      Lucee::Region<CDIM, int> localRgn = momentGlobal.getRegion();
      Lucee::Region<CDIM, int> localExtRgn = momentGlobal.getExtRegion();
      
      int lowerConf[CDIM], upperConf[CDIM], lg[CDIM], ug[CDIM];
      for (int i=0; i<CDIM; ++i)
      {
        lowerConf[i] = localRgn.getLower(i);
        upperConf[i] = localRgn.getUpper(i);
        lg[i] = localRgn.getLower(i) - localExtRgn.getLower(i);
        ug[i] = localExtRgn.getUpper(i) - localRgn.getUpper(i);
      }
      Lucee::Region<CDIM, int> rgnConf(lowerConf, upperConf);
      momentLocal = new Lucee::Field<CDIM, double>(rgnConf, momentGlobal.getNumComponents(), lg, ug);
    }
    
    // local phase-space region to update
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    int idx[NDIM];
    double xc[PNC];
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    unsigned nlocalConf = confBasis->getNumNodes();
    unsigned nlocalPhase = phaseBasis->getNumNodes();

    // total number of local configuration space cells
    int localPositionCells = momentGlobal.getExtRegion().getVolume();

    // clear out contents of output field
    (*momentLocal) = 0.0;

    // iterators into fields
    Lucee::ConstFieldPtr<double> distFPtr = distF.createConstPtr();
    Lucee::FieldPtr<double> momentPtr = momentLocal->createPtr();
    Eigen::VectorXd distfVec(nlocalPhase);

    std::vector<Eigen::VectorXd> resultVector = std::vector<Eigen::VectorXd>(nMom);
    for (int i=0; i<nMom; ++i)
      resultVector[i] = Eigen::VectorXd::Zero(nlocalConf);    
   
    // Temporary array to make computation of particle energy cleaner
    std::vector<Eigen::VectorXd> resultVectorEnergy = std::vector<Eigen::VectorXd>(VDIM);
    for (int i=0; i<VDIM; ++i)
      resultVectorEnergy[i] = Eigen::VectorXd::Zero(nlocalConf);  

    // accumulate contribution to moment from this cell
    while(seq.step())
    {
      seq.fillWithIndex(idx);
      grid.setIndex(idx);
      grid.getCentroid(xc);
      momentLocal->setPtr(momentPtr, idx);
      distF.setPtr(distFPtr, idx);

      for (int i=0; i<nlocalPhase; ++i)
        distfVec(i) = distFPtr[i];

      if (calcMom == 0)
      {
        // density
        resultVector[0].noalias() = mom0Matrix*distfVec;
      }
      else if (calcMom == 1)
      { // momentum
        for (int h=0; h<VDIM; ++h)
        {
          int momDir = h+CDIM;
          resultVector[h].noalias() = (mom1Matrix[h] + xc[momDir]*mom0Matrix)*distfVec;
        }
      }
      else if (calcMom == 2)
      {
          resultVector[0] *= 0.0;
          int ctr = 0;
          for (int h=0; h<VDIM; ++h)
          {
            int momDir = h+CDIM;
            resultVectorEnergy[h].noalias() = (mom2Matrix[ctr]+ 2*xc[momDir]*mom1Matrix[h] + xc[momDir]*xc[momDir]*mom0Matrix)*distfVec;
            ctr = ctr+VDIM-h;
            resultVector[0] += resultVectorEnergy[h];
          }
          resultVector[0] *= 0.5;
      }

      for (int i=0; i<nlocalConf; ++i)
        for (int j=0; j<nMom; ++j)
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
    momComm->allreduce(xsize, &momentLocal->first(), &momentGlobal.first(), TX_SUM);

    return Lucee::UpdaterStatus();
  }

  template <unsigned CDIM, unsigned VDIM>  
  void
  BoltzmannCollisionOperator<CDIM, VDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<CDIM+VDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<CDIM, double>));
  }

  template <unsigned CDIM, unsigned VDIM>
  void
  BoltzmannCollisionOperator<CDIM, VDIM>::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); ++rowIndex)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); ++colIndex)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }

// instantiations
  template class BoltzmannCollisionOperator<1,1>;
  template class BoltzmannCollisionOperator<1,2>;
  template class BoltzmannCollisionOperator<1,3>;
  template class BoltzmannCollisionOperator<2,2>;
  template class BoltzmannCollisionOperator<2,3>;  
}
