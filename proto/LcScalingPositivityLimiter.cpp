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

    // flag for storing energy
    storeEnergy = false;
    if (tbl.hasBool("storeEnergy"))
      storeEnergy = tbl.getBool("storeEnergy");
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

    const Lucee::Field<NDIM, double>& distfIn =
      this->getInp<Lucee::Field<NDIM, double> >(0);
    Lucee::Field<NDIM, double>&       distfOut =
      this->getOut<Lucee::Field<NDIM, double> >(0);

    Lucee::Region<NDIM, int>    localRgn  = grid.getLocalRegion();
    Lucee::Region<NDIM, double> compSpace = grid.getComputationalSpace();

    Lucee::ConstFieldPtr<double> distfInPtr  = distfIn.createConstPtr();
    Lucee::FieldPtr<double>      distfOutPtr = distfOut.createPtr();

    Lucee::Matrix<double> phaseNodeCoords(phaseBasis->getNumNodes(), PNC);

    distfOut = 0.0; // use distfOut to store increment initially    
    
    int idx[NDIM];
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);

    unsigned numNodesPhase = phaseBasis->getNumNodes();

    std::vector<double> weights(numNodesPhase);
    phaseBasis->getWeights(weights);
    
    double cellAverage, cellMin, invCellVolume;
    double theta;

    invCellVolume = 1.0/pow(2, NDIM);
    
    while (seq.step()) {
      seq.fillWithIndex(idx);
      distfIn.setPtr(distfInPtr, idx);
      distfOut.setPtr(distfInPtr, idx);
      
      phaseBasis->setIndex(idx);
      
      cellMin = 1e9;
      cellAverage = 0;
      for (unsigned nodeIdx = 0; nodeIdx<numNodesPhase; ++nodeIdx) {
	cellAverage += weights[nodeIdx]*distfInPtr[nodeIdx];
	if (distfInPtr[nodeIdx] < cellMin) {
	  cellMin = distfInPtr[nodeIdx];
	}
      }
      cellAverage *= invCellVolume;
      theta = std::min(1.0, cellAverage/(cellAverage-cellMin));

      for (unsigned nodeIdx = 0; nodeIdx<numNodesPhase; ++nodeIdx) {
	distfOutPtr[nodeIdx] = cellAverage + 
	  theta*(distfInPtr[nodeIdx] - cellAverage);
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

