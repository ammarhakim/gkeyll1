/**
 * @file	LcNodalContinuumKineticSEE
 *
 * @brief	Updater for calculation Secondary Emitted Electron (SEE) BC
 */

// lucee includes
#include <LcNodalContinuumKineticSEE.h>

namespace Lucee
{
  static const unsigned LC_LOWER_EDGE = 0;
  static const unsigned LC_UPPER_EDGE = 1;

  // set ids for module system
  template <> const char *NodalContinuumKineticSEE<1, 1>::id =
    "NodalContinuumKineticSEE1X1V";
  template <> const char *NodalContinuumKineticSEE<1, 2>::id =
    "NodalContinuumKineticSEE1X2V";
  template <> const char *NodalContinuumKineticSEE<1, 3>::id =
    "NodalContinuumKineticSEE1X3V";
  template <> const char *NodalContinuumKineticSEE<2, 2>::id =
    "NodalContinuumKineticSEE2X2V";
  template <> const char *NodalContinuumKineticSEE<2, 3>::id =
    "NodalContinuumKineticSEE2X3V";
  //template <> const char *NodalContinuumKineticSEE<3, 3>::id = "NodalContinuumKineticSEE3X3V";

  template <unsigned CDIM, unsigned VDIM>
  bool
  NodalContinuumKineticSEE<CDIM,VDIM>::sameConfigCoords(unsigned n, unsigned cn, double dxMin,
    const Lucee::Matrix<double>& phaseC, const Lucee::Matrix<double>& confC)
  {
    for (unsigned dim = 0; dim<CDIM; ++dim)
      if (! (std::fabs(phaseC(n, dim)-confC(cn, dim))<1e-4*dxMin) )
        return false;
    return true;
  }

  //------------------------------------------------------------------
  //-- Constructor and Destructor ------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  NodalContinuumKineticSEE<CDIM,VDIM>::NodalContinuumKineticSEE()
    : UpdaterIfc()
  {
  }
  template <unsigned CDIM, unsigned VDIM>
  NodalContinuumKineticSEE<CDIM, VDIM>::~NodalContinuumKineticSEE()
  {
  }

  //------------------------------------------------------------------
  //-- Read input ----------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  void
  NodalContinuumKineticSEE<CDIM, VDIM>::readInput(Lucee::LuaTable& tbl)
  {
    const unsigned NDIM = CDIM+VDIM;
    // call base class method
    UpdaterIfc::readInput(tbl);

    // get hold on the basis
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("phaseBasis"))
      phaseBasis = 
	&tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("phaseBasis");
    else
      throw Lucee::Except("NodalContinuumKineticSEE::readInput: Must specify phase space basis to use using 'phaseBasis'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<CDIM> >("confBasis"))
      confBasis =
	&tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<CDIM> >("confBasis");
    else
      throw Lucee::Except("NodalContinuumKineticSEE::readInput: Must specify configuration space basis to use using 'confBasis'");

    if (tbl.hasString("edge")) {
      std::string edgeStr = tbl.getString("edge");
      if (edgeStr == "lower")
	edge = LC_LOWER_EDGE;
      else
	edge = LC_UPPER_EDGE;
    } else
      throw Lucee::Except("NodalContinuumKineticSEE::readInput: Must specify edge using 'edge'");

    if (CDIM > 1) {
      if (tbl.hasNum("dir"))
	dir = &tbl.getNum("dir");
      else
	throw Lucee::Except("NodalContinuumKineticSEE::readInput: Must specify direction using 'dir'");
    } else {
      dir = 0;
    }
  }

  //------------------------------------------------------------------
  //-- Initialize ----------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  void 
  NodalContinuumKineticSEE<CDIM, VDIM>::initialize()
  {
    UpdaterIfc::initialize();
 
    const unsigned NDIM = CDIM+VDIM;
    /*unsigned nlocal = phaseBasis->getNumNodes();

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
          "NodalContinuumKineticSEE::initialize: No matching configuration space node for phase-space node ");
        lce << n;
        throw lce;
      }
      }*/
  }

  //------------------------------------------------------------------
  //-- Update --------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  Lucee::UpdaterStatus 
  NodalContinuumKineticSEE<CDIM, VDIM>::update(double t)
  {
    const unsigned NDIM = CDIM+VDIM;

    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    Lucee::Field<NDIM, double>& f =
      this->getOut<Lucee::Field<NDIM, double> >(0);
    Lucee::FieldPtr<double> fPtrEdge = f.createPtr();
    Lucee::FieldPtr<double> fPtrGhost = f.createPtr();

    Lucee::Matrix<double> phaseNodeCoords(phaseBasis->getNumNodes(), PNC);  

    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    Lucee::Region<NDIM, int> externalRgn = f.getExtRegion();

    int lowerConf[CDIM], upperConf[CDIM];
    for (int i = 0; i < VDIM; i++) {
      lowerVel[i] = localRgn.getLower(CDIM+i);
      upperVel[i] = localRgn.getUpper(CDIM+i);
    }
    Lucee::Region<VDIM, int> velRgn(lower, upper);
    int lowerVel[VDIM], upperVel[VDIM];
    for (int i = 0; i < VDIM; i++) {
      lowerVel[i] = localRgn.getLower(CDIM+i);
      upperVel[i] = localRgn.getUpper(CDIM+i);
    }
    Lucee::Region<VDIM, int> velRgn(lower, upper);
    
    int idx[NDIM], idxIn[VDIM], idxOut[VDIM];
    Lucee::RowMajorSequencer<NDIM> seqVelIn(velRgn);
    Lucee::RowMajorSequencer<NDIM> seqVelOut(velRgn);

    unsigned numNodesPhase = phaseBasis->getNumNodes();
    
    while (seqOut.step()) {
      seq.fillWithIndex(idx);
      q.setPtr(qPtr, idx);
      zerothMoment.setPtr(zerothMomentPtr, idx);
      firstMoment.setPtr(firstMomentPtr, idx);
      secondMoment.setPtr(secondMomentPtr, idx);
      
      phaseBasis->setIndex(idx);
      phaseBasis->getNodalCoordinates(phaseNodeCoords);
      
      for (unsigned nodeIdx = 0; nodeIdx<numNodesConf; ++nodeIdx) {
	dens[nodeIdx] = zerothMomentPtr[nodeIdx];
	invDens[nodeIdx] = 1/dens[nodeIdx];
	for (unsigned dim = 0; dim<VDIM; ++dim) {
	  vDrift(nodeIdx, dim) = invDens[nodeIdx] *
	    firstMomentPtr[nodeIdx*VDIM + dim];
	  
	  double temp = invDens[nodeIdx] *
	    secondMomentPtr[nodeIdx*VDIM + dim] -
	    vDrift(nodeIdx, dim)*vDrift(nodeIdx, dim);
	  if (temp < 0)
	    return Lucee::UpdaterStatus(false, 0);
	  invVTerm2(nodeIdx, dim) = 1/temp;
	  invVTerm(nodeIdx, dim) = sqrt(invVTerm2(nodeIdx, dim));
	}  	
      }
      for (unsigned nodeIdx = 0; nodeIdx<numNodesPhase; ++nodeIdx) {
	// get number density
	if (!arbitraryDensity)
	  n = dens[phaseConfMap[nodeIdx]];
	else {
	  const Lucee::Field<CDIM, double>& zerothMomentOut =
	    this->getInp<Lucee::Field<CDIM, double> >(3);
	  Lucee::ConstFieldPtr<double> zerothMomentOutPtr = 
	    zerothMomentOut.createConstPtr();
	  zerothMomentOut.setPtr(zerothMomentOutPtr, idx);
	  n = zerothMomentOutPtr[phaseConfMap[nodeIdx]];	  
	}
	
	// get volocity for Maxwellian
	if (zeroDriftOutput) 
	  for (unsigned dim = 0; dim<VDIM; ++dim)
	    v[dim] = phaseNodeCoords(nodeIdx, CDIM+dim);
	else
	  for (unsigned dim = 0; dim<VDIM; ++dim)
	    v[dim] = phaseNodeCoords(nodeIdx, CDIM+dim)-
	      vDrift(phaseConfMap[nodeIdx], dim);
	
	// get thermal velocity
	for (unsigned dim = 0; dim<VDIM; ++dim) {
	  invVt2[dim] = invVTerm2(phaseConfMap[nodeIdx], dim);
	  invVt[dim] = invVTerm(phaseConfMap[nodeIdx], dim);
	}
	
	qPtr[nodeIdx] = evaluateMaxwell(n, v, invVt, invVt2);
      }
    }
    
    return Lucee::UpdaterStatus(true, 0);
  }

  //------------------------------------------------------------------
  //-- Declare types -------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  void NodalContinuumKineticSEE<CDIM, VDIM>::declareTypes()
  {
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }
 
  //------------------------------------------------------------------
  //-- Instantiations ------------------------------------------------
  template class NodalContinuumKineticSEE<1, 1>;
  template class NodalContinuumKineticSEE<1, 2>;
  template class NodalContinuumKineticSEE<1, 3>;
  template class NodalContinuumKineticSEE<2, 2>;
  template class NodalContinuumKineticSEE<2, 3>;
  template class NodalContinuumKineticSEE<3, 3>;
}

