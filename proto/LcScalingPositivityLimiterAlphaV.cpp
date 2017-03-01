/**
 * @file	LcScalingPositivityLimiterAlphaV.cpp
 *
 * @brief	Calculate the drag factor to remove energy added by the 
 *              ScalingpositivityLimiter
 */

// lucee includes
#include <LcScalingPositivityLimiterAlphaV.h>

namespace Lucee
{
// set ids for module system
  template <> const char *ScalingPositivityLimiterAlphaV<1, 1>::id = 
    "ScalingPositivityLimiterAlphaV1X1V";
  template <> const char *ScalingPositivityLimiterAlphaV<1, 2>::id =
    "ScalingPositivityLimiterAlphaV1X2V";
  template <> const char *ScalingPositivityLimiterAlphaV<1, 3>::id = 
    "ScalingPositivityLimiterAlphaV1X3V";
  template <> const char *ScalingPositivityLimiterAlphaV<2, 2>::id =
    "ScalingPositivityLimiterAlphaV2X2V";
  template <> const char *ScalingPositivityLimiterAlphaV<2, 3>::id =
    "ScalingPositivityLimiterAlphaV2X3V";
  /*template <> const char *ScalingPositivityLimiterAlphaV<3, 3>::id =
    "ScalingPositivityLimiterAlphaV3X3V";*/

  template <unsigned CDIM, unsigned VDIM>
  bool
  ScalingPositivityLimiterAlphaV<CDIM,VDIM>::sameConfigCoords(unsigned n, unsigned cn, double dxMin,
    const Lucee::Matrix<double>& phaseC, const Lucee::Matrix<double>& confC)
  {
    for (unsigned dim = 0; dim<CDIM; ++dim)
      if (! (std::fabs(phaseC(n, dim)-confC(cn, dim))<1e-4*dxMin) )
        return false;
    return true;
  }
  //----------------------------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  ScalingPositivityLimiterAlphaV<CDIM,VDIM>::ScalingPositivityLimiterAlphaV()
    : UpdaterIfc()
  {
  }
  template <unsigned CDIM, unsigned VDIM>
  ScalingPositivityLimiterAlphaV<CDIM, VDIM>::~ScalingPositivityLimiterAlphaV()
  {
  }
  //----------------------------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  void
  ScalingPositivityLimiterAlphaV<CDIM, VDIM>::readInput(Lucee::LuaTable& tbl)
  {
    const unsigned NDIM = CDIM+VDIM;
    // call base class method
    UpdaterIfc::readInput(tbl);

    // get hold on the basis
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("phaseBasis"))
      phaseBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("phaseBasis");
    else
      throw Lucee::Except("ScalingPositivityLimiterAlphaV::readInput: Must specify phase space basis to use using 'phaseBasis'");
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<CDIM> >("confBasis"))
      confBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<CDIM> >("confBasis");
    else
      throw Lucee::Except("ScalingPositivityLimiterAlphaV::readInput: Must specify configuration space basis to use using 'confBasis'");
  }
  //----------------------------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  void 
  ScalingPositivityLimiterAlphaV<CDIM, VDIM>::initialize()
  {
    UpdaterIfc::initialize();
 
    const unsigned NDIM = CDIM+VDIM;
    //unsigned nlocal = phaseBasis->getNumNodes();

    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[NDIM];
    seq.fillWithIndex(idx);
    phaseBasis->setIndex(idx);
    confBasis->setIndex(idx); // only first CDIM elements are used

    // compute mapping of phase-space nodes to configuration space
    // nodes. The assumption here is that the node layout in phase-space
    // and configuration space are such that each node in phase-space has
    // exactly one node co-located with it in configuration space. No
    // "orphan" phase-space node are allowed, and an exception is thrown
    // if that occurs.
    /*
    phaseConfMap.resize(nlocal);
    Lucee::Matrix<double> phaseNodeCoords(phaseBasis->getNumNodes(), PNC);
    Lucee::Matrix<double> confNodeCoords(confBasis->getNumNodes(), CNC);

    double dxMin = grid.getDx(0);
    for (unsigned d=1; d<CDIM; ++d)
      dxMin = std::min(dxMin, grid.getDx(d));

    phaseBasis->getNodalCoordinates(phaseNodeCoords);
    confBasis->getNodalCoordinates(confNodeCoords);
    for (unsigned n=0; n<nlocal; ++n)
    {
      bool pcFound = false;
      for (unsigned cn=0; cn<confBasis->getNumNodes(); ++cn)
        if (sameConfigCoords(n, cn, dxMin, phaseNodeCoords, confNodeCoords))
        {
          phaseConfMap[n] = cn;
          pcFound = true;
          break;
        }
      if (!pcFound)
      {
        Lucee::Except lce(
          "ScalingPositivityLimiterAlphaV::initialize: No matching configuration space node for phase-space node ");
        lce << n;
        throw lce;
      }
      }*/
  }
  //----------------------------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  Lucee::UpdaterStatus 
  ScalingPositivityLimiterAlphaV<CDIM, VDIM>::update(double t)
  {
    const unsigned NDIM = CDIM+VDIM;

    double dt = t-this->getCurrTime();

    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    const Lucee::Field<CDIM, double>& energyOld =
      this->getInp<Lucee::Field<CDIM, double> >(0);
    const Lucee::Field<CDIM, double>& energyNew =
      this->getInp<Lucee::Field<CDIM, double> >(1);
      
    Lucee::Field<NDIM, double>& alphaV =
      this->getOut<Lucee::Field<NDIM, double> >(0);

    Lucee::Region<NDIM, int>    localRgn  = grid.getLocalRegion();
    Lucee::Region<NDIM, double> compSpace = grid.getComputationalSpace();

    Lucee::ConstFieldPtr<double> energyOldPtr =
      energyOld.createConstPtr();
    Lucee::ConstFieldPtr<double> energyNewPtr =
      energyNew.createConstPtr();

    Lucee::FieldPtr<double> alphaVPtr = alphaV.createPtr();

    Lucee::Matrix<double> phaseNodeCoords(phaseBasis->getNumNodes(), PNC);

    alphaV = 0.0; // use q to store increment initially    
    
    int idx[NDIM];
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);

    unsigned numNodesPhase = phaseBasis->getNumNodes();
    unsigned numNodesConf  = confBasis->getNumNodes();

    std::vector<double> weights(numNodesConf);
    confBasis->getWeights(weights);

    double cellVolume, invCellVolume = 0;
    for (unsigned nodeIdx = 0; nodeIdx<numNodesConf; ++nodeIdx)
      cellVolume += weights[nodeIdx];
    invCellVolume = 1.0/cellVolume;

    double alpha;
    
    while (seq.step()) {
      seq.fillWithIndex(idx);
      energyOld.setPtr(energyOldPtr, idx);
      energyNew.setPtr(energyNewPtr, idx);
      alphaV.setPtr(alphaVPtr, idx);
      
      phaseBasis->setIndex(idx);
      phaseBasis->getNodalCoordinates(phaseNodeCoords);
      
      alpha = 0;
      for (unsigned nodeIdx = 0; nodeIdx<numNodesConf; ++nodeIdx)
	alpha += weights[nodeIdx]*0.5*
	  (energyNewPtr[nodeIdx] - energyOldPtr[nodeIdx])/
	  (dt*energyOldPtr[nodeIdx]);
      alpha *= invCellVolume;

      for (unsigned nodeIdx = 0; nodeIdx<numNodesPhase; ++nodeIdx) 
	for (unsigned dim = 0; dim<VDIM; ++dim)
	  alphaVPtr[nodeIdx*VDIM + dim] = 
	    alpha*phaseNodeCoords(nodeIdx, CDIM+dim);
    }
    
    return Lucee::UpdaterStatus(true, 0);
  }
  //----------------------------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  void ScalingPositivityLimiterAlphaV<CDIM, VDIM>::declareTypes()
  {
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }
  //----------------------------------------------------------------------------
  // instantiations
  template class ScalingPositivityLimiterAlphaV<1, 1>;
  template class ScalingPositivityLimiterAlphaV<1, 2>;
  template class ScalingPositivityLimiterAlphaV<1, 3>;
  template class ScalingPositivityLimiterAlphaV<2, 2>;
  template class ScalingPositivityLimiterAlphaV<2, 3>;
  //template class ScalingPositivityLimiterAlphaV<3, 3>;
}

