/**
 * @file	LcScalingPositivityLimiter.cpp
 *
 * @brief	A scaling limiter that ensures positivity of nodal values
 */

// lucee includes
#include <LcScalingPositivityLimiter.h>

namespace Lucee
{
// set ids for module system
  template <> const char *ScalingPositivityLimiter<1, 1>::id =
    "ScalingPositivityLimiter1X1V";
  template <> const char *ScalingPositivityLimiter<1, 2>::id =
    "ScalingPositivityLimiter1X2V";
  template <> const char *ScalingPositivityLimiter<1, 3>::id = 
    "ScalingPositivityLimiter1X3V";
  template <> const char *ScalingPositivityLimiter<2, 2>::id =
    "ScalingPositivityLimiter2X2V";
  template <> const char *ScalingPositivityLimiter<2, 3>::id = 
    "ScalingPositivityLimiter2X3V";
  //template <> const char *ScalingPositivityLimiter<3, 3>::id = "ScalingPositivityLimiter3X3V";

  //----------------------------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  ScalingPositivityLimiter<CDIM,VDIM>::ScalingPositivityLimiter()
    : UpdaterIfc()
  {
  }
  template <unsigned CDIM, unsigned VDIM>
  ScalingPositivityLimiter<CDIM, VDIM>::~ScalingPositivityLimiter()
  {
  }
  //----------------------------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  void
  ScalingPositivityLimiter<CDIM, VDIM>::readInput(Lucee::LuaTable& tbl)
  {
    const unsigned NDIM = CDIM+VDIM;
    // call base class method
    UpdaterIfc::readInput(tbl);

    // get hold on the basis
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      phaseBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("ScalingPositivityLimiter::readInput: Must specify phase space basis to use using 'basis'");
  }
  //----------------------------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  void 
  ScalingPositivityLimiter<CDIM, VDIM>::initialize()
  {
    UpdaterIfc::initialize();
 
    const unsigned NDIM = CDIM+VDIM;
    unsigned nlocal = phaseBasis->getNumNodes();

    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[NDIM];
    seq.fillWithIndex(idx);
    phaseBasis->setIndex(idx);
  }
  //----------------------------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  Lucee::UpdaterStatus 
  ScalingPositivityLimiter<CDIM, VDIM>::update(double t)
  {
    const unsigned NDIM = CDIM+VDIM;

    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    Lucee::Field<NDIM, double>& distf =
      this->getOut<Lucee::Field<NDIM, double> >(0);
    Lucee::FieldPtr<double> distfPtr = distf.createPtr();

    Lucee::Region<NDIM, int>    localRgn  = grid.getLocalRegion();
    Lucee::Region<NDIM, double> compSpace = grid.getComputationalSpace();
    
    int idx[NDIM];
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);

    unsigned numNodesPhase = phaseBasis->getNumNodes();

    std::vector<double> weights(numNodesPhase);
    phaseBasis->getWeights(weights);
    
    double cellAvg, cellMin, invCellVolume;
    double theta;

    invCellVolume = 1.0/grid.getVolume();
    
    while (seq.step()) {
      seq.fillWithIndex(idx);
      distf.setPtr(distfPtr, idx);
      
      phaseBasis->setIndex(idx);
            
      cellMin = 1e9;
      cellAvg = 0;
      for (unsigned nodeIdx = 0; nodeIdx<numNodesPhase; ++nodeIdx) {
	cellAvg += weights[nodeIdx]*distfPtr[nodeIdx];
	if (distfPtr[nodeIdx] < cellMin) {
	  cellMin = distfPtr[nodeIdx];
	}
      }
      cellAvg *= invCellVolume;
      theta = fmin(1.0, cellAvg/(cellAvg-cellMin));

      for (unsigned nodeIdx = 0; nodeIdx<numNodesPhase; ++nodeIdx) {
	distfPtr[nodeIdx] = cellAvg + 
	  theta*(distfPtr[nodeIdx] - cellAvg);
      }		  
    }
    return Lucee::UpdaterStatus(true, 0);
  }
  //----------------------------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  void ScalingPositivityLimiter<CDIM, VDIM>::declareTypes()
  {
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }
  //----------------------------------------------------------------------------
  // instantiations
  template class ScalingPositivityLimiter<1, 1>;
  template class ScalingPositivityLimiter<1, 2>;
  template class ScalingPositivityLimiter<1, 3>;
  template class ScalingPositivityLimiter<2, 2>;
  template class ScalingPositivityLimiter<2, 3>;
  //template class ScalingPositivityLimiter<3, 3>;
}

