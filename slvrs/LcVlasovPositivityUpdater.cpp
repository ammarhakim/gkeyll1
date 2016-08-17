/**
 * @file	LcVlasovPositivityUpdater.cpp
 *
 * @brief	Updater to enforce positivity preservation of the distribution function in Vlasov simulations
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAlignedRectCoordSys.h>
#include <LcField.h>
#include <LcLinAlgebra.h>
#include <LcMathLib.h>
#include <LcVlasovPositivityUpdater.h>
#include <LcStructuredGridBase.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
// set id for module system
  template <> const char *VlasovPositivityUpdater<1,1>::id = "VlasovPositivity1X1V";
  template <> const char *VlasovPositivityUpdater<1,2>::id = "VlasovPositivity1X2V";
  template <> const char *VlasovPositivityUpdater<1,3>::id = "VlasovPositivity1X3V";
  template <> const char *VlasovPositivityUpdater<2,2>::id = "VlasovPositivity2X2V";
  template <> const char *VlasovPositivityUpdater<2,3>::id = "VlasovPositivity2X3V";
  //template <> const char *VlasovPositivityUpdater<3,3>::id = "VlasovPositivity3X3V";


  template <unsigned CDIM, unsigned VDIM>
  bool
  VlasovPositivityUpdater<CDIM,VDIM>::sameConfigCoords(unsigned n, unsigned cn, double dxMin,
    const Eigen::MatrixXd& phaseC, const Eigen::MatrixXd& confC)
  {
    for (unsigned dim = 0; dim<CDIM; ++dim)
      if (! (std::fabs(phaseC(n, dim)-confC(cn, dim))<1e-4*dxMin) )
        return false;
    return true;
  }

  template <unsigned CDIM, unsigned VDIM>
  VlasovPositivityUpdater<CDIM, VDIM>::VlasovPositivityUpdater()
    : Lucee::UpdaterIfc()
  {
  }

  template <unsigned CDIM, unsigned VDIM>  
  void
  VlasovPositivityUpdater<CDIM, VDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);    
    const unsigned NDIM = CDIM+VDIM;

    // get hold of phase space element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("phaseBasis"))
      phaseBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("phaseBasis");
    else
      throw Lucee::Except(
        "VlasovPositivityUpdater::readInput: Must specify phase-space basis using 'phaseBasis'");

    // get hold of configuration space element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<CDIM> >("confBasis"))
      confBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<CDIM> >("confBasis");
    else
      throw Lucee::Except(
        "VlasovPositivityUpdater::readInput: Must specify configuration-space basis using 'confBasis'");
  }

  template <unsigned CDIM, unsigned VDIM>
  void
  VlasovPositivityUpdater<CDIM, VDIM>::initialize()
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
          "VlasovPositivityUpdater::initialize: No matching configuration space quadrature node for phase-space quadrature node ");
        lce << n;
        throw lce;
      }
    }

    // Get configuration space mass matrix
    Lucee::Matrix<double> tempMassMatrixConf(nlocalConf, nlocalConf);
    confBasis->getMassMatrix(tempMassMatrixConf);
    Eigen::MatrixXd massMatrixConf(nlocalConf, nlocalConf);
    copyLuceeToEigen(tempMassMatrixConf, massMatrixConf);

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
  }
  template <unsigned CDIM, unsigned VDIM>
  Lucee::UpdaterStatus
  VlasovPositivityUpdater<CDIM, VDIM>::update(double t)
  {
    const unsigned NDIM = CDIM+VDIM;
    // get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    // get output field (CDIM + VDIM) 
    Lucee::Field<NDIM, double>& distfOut = this->getOut<Lucee::Field<NDIM, double> >(0);

    unsigned nlocalConf = confBasis->getNumNodes();
    unsigned nlocalPhase = phaseBasis->getNumNodes();
    Lucee::FieldPtr<double> distfPtr = distfOut.createPtr();

    int idx[NDIM];
    // local region to update
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);

    // Used in the sequencer
    Eigen::VectorXd distfVector(nlocalPhase);

    while (seq.step())
    {
      seq.fillWithIndex(idx);
      grid.setIndex(idx);
      distfOut.setPtr(distfPtr, idx);

      // Fill out distfVector
      for (int i=0; i<nlocalPhase; ++i)
        distfVector(i) = distfPtr[i];

      // Compute density
      double originalNum = (mom0Matrix*distfVector).sum();
      if (originalNum<0.0)
      {
        //std::cout << "entire cell negative (density = " << originalNum << ")" << std::endl;
        for (int i=0; i<nlocalPhase; ++i)
        {
          if (distfPtr[i] < 0.0)
            distfPtr[i] = 0.0;
        }
        continue;
      }
      else if (originalNum==0.0)
      {
        //std::cout << "(positivity) cell is zero. skipping." << std::endl;
        continue;
      }
      // Zero out distfVector entries that are negative
      for (int i=0; i<nlocalPhase; ++i)
      {
        if (distfVector(i)<0.0)
          distfVector(i) = 0.0;
      }

      // Compute modified density. This will be greater than originalNum
      double modifiedNum = (mom0Matrix*distfVector).sum();

      if (modifiedNum==0.0)
      {
        //std::cout << "New cell is all zero" << std::endl;
        for (int i=0; i<nlocalPhase; ++i)
          distfPtr[i] = 0.0;
      }
      else
      {
        // Write modified values to distfOut
        for (int i=0; i<nlocalPhase; ++i)
          distfPtr[i] = (originalNum/modifiedNum)*distfVector(i);
      }
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned CDIM, unsigned VDIM>  
  void
  VlasovPositivityUpdater<CDIM, VDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<CDIM+VDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<CDIM, double>));
  }

  template <unsigned CDIM, unsigned VDIM>
  void
  VlasovPositivityUpdater<CDIM, VDIM>::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); ++rowIndex)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); ++colIndex)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }

  // instantiations
  template class VlasovPositivityUpdater<1,1>;
  template class VlasovPositivityUpdater<1,2>;
  template class VlasovPositivityUpdater<1,3>;
  template class VlasovPositivityUpdater<2,2>;
  template class VlasovPositivityUpdater<2,3>;
}
