/**
 * @file	LcNodalDgScalingLimiterUpdater.cpp
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
#include <LcNodalDgScalingLimiterUpdater.h>
#include <LcStructuredGridBase.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
// set id for module system
  template <> const char *NodalDgScalingLimiterUpdater<1,1>::id = "NodalDgScalingLimiter1X1V";
  template <> const char *NodalDgScalingLimiterUpdater<1,2>::id = "NodalDgScalingLimiter1X2V";
  template <> const char *NodalDgScalingLimiterUpdater<1,3>::id = "NodalDgScalingLimiter1X3V";
  template <> const char *NodalDgScalingLimiterUpdater<2,2>::id = "NodalDgScalingLimiter2X2V";
  template <> const char *NodalDgScalingLimiterUpdater<2,3>::id = "NodalDgScalingLimiter2X3V";
  //template <> const char *NodalDgScalingLimiterUpdater<3,3>::id = "NodalDgScalingLimiter3X3V";

  template <unsigned CDIM, unsigned VDIM>
  NodalDgScalingLimiterUpdater<CDIM, VDIM>::NodalDgScalingLimiterUpdater()
    : Lucee::UpdaterIfc()
  {
  }

  template <unsigned CDIM, unsigned VDIM>  
  void
  NodalDgScalingLimiterUpdater<CDIM, VDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);    
    const unsigned NDIM = CDIM+VDIM;

    // get hold of phase space element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("phaseBasis"))
      phaseBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("phaseBasis");
    else
      throw Lucee::Except(
        "NodalDgScalingLimiterUpdater::readInput: Must specify phase-space basis using 'phaseBasis'");
  }

  template <unsigned CDIM, unsigned VDIM>
  void
  NodalDgScalingLimiterUpdater<CDIM, VDIM>::initialize()
  {
    Lucee::UpdaterIfc::initialize();    
    const unsigned NDIM = CDIM+VDIM;

    // get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    // local region to update
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    // get number of nodes in configuration space and phase space
    unsigned nlocalPhase = phaseBasis->getNumNodes();
    // get weights for quadrature
    volWeights = Eigen::VectorXd::Zero(nlocalPhase);
    std::vector<double> weights(nlocalPhase);
    phaseBasis->getWeights(weights);

    double vol = grid.getVolume();
    // normalize weights as we are computing averages and not integrals over cell
    for (unsigned i=0; i<nlocalPhase; ++i)
      volWeights(i) = weights[i]/vol;
  }
  template <unsigned CDIM, unsigned VDIM>
  Lucee::UpdaterStatus
  NodalDgScalingLimiterUpdater<CDIM, VDIM>::update(double t)
  {
    const unsigned NDIM = CDIM+VDIM;
    // get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    // get output field (CDIM + VDIM) 
    Lucee::Field<NDIM, double>& distfOut = this->getOut<Lucee::Field<NDIM, double> >(0);

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

      // Compute number of particles
      double originalNum = (distfVector.cwiseProduct(volWeights)).sum();
      if (originalNum<0.0)
      {
        std::cout << "entire cell negative (density = " << originalNum << ")" << std::endl;
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
      double modifiedNum = (distfVector.cwiseProduct(volWeights)).sum();

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
  NodalDgScalingLimiterUpdater<CDIM, VDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<CDIM+VDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<CDIM, double>));
  }

  template <unsigned CDIM, unsigned VDIM>
  void
  NodalDgScalingLimiterUpdater<CDIM, VDIM>::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); ++rowIndex)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); ++colIndex)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }

  // instantiations
  template class NodalDgScalingLimiterUpdater<1,1>;
  template class NodalDgScalingLimiterUpdater<1,2>;
  template class NodalDgScalingLimiterUpdater<1,3>;
  template class NodalDgScalingLimiterUpdater<2,2>;
  template class NodalDgScalingLimiterUpdater<2,3>;
}
