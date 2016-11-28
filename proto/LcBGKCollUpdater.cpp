/**
 * @file	LcBGKCollUpdater.cpp
 *
 * @brief	Updater for applying the BGK update from moments
 */

// lucee includes
#include <LcBGKCollUpdater.h>

namespace Lucee
{
// set ids for module system
  template <> const char *BGKCollUpdater<1, 1>::id = "BGKCollUpdater1X1V";
  template <> const char *BGKCollUpdater<1, 2>::id = "BGKCollUpdater1X2V";
  template <> const char *BGKCollUpdater<1, 3>::id = "BGKCollUpdater1X3V";
  template <> const char *BGKCollUpdater<2, 2>::id = "BGKCollUpdater2X2V";
  template <> const char *BGKCollUpdater<2, 3>::id = "BGKCollUpdater2X3V";
  //template <> const char *BGKCollUpdater<3, 3>::id = "BGKCollUpdater3X3V";

  template <unsigned CDIM, unsigned VDIM>
  bool
  BGKCollUpdater<CDIM,VDIM>::sameConfigCoords(unsigned n, unsigned cn, double dxMin,
    const Lucee::Matrix<double>& phaseC, const Lucee::Matrix<double>& confC)
  {
    for (unsigned dim = 0; dim<CDIM; ++dim)
      if (! (std::fabs(phaseC(n, dim)-confC(cn, dim))<1e-4*dxMin) )
        return false;
    return true;
  }
  //----------------------------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  BGKCollUpdater<CDIM,VDIM>::BGKCollUpdater()
    : UpdaterIfc()
  {
    maxwellNorm  = 1/sqrt(2*M_PI);
  }
  template <unsigned CDIM, unsigned VDIM>
  BGKCollUpdater<CDIM, VDIM>::~BGKCollUpdater()
  {
  }
  //----------------------------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  void
  BGKCollUpdater<CDIM, VDIM>::readInput(Lucee::LuaTable& tbl)
  {
    const unsigned NDIM = CDIM+VDIM;
    // call base class method
    UpdaterIfc::readInput(tbl);

    // get hold on the basis
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("phaseBasis"))
      phaseBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("phaseBasis");
    else
      throw Lucee::Except("BGKCollUpdater::readInput: Must specify phase space basis to use using 'phaseBasis'");
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<CDIM> >("confBasis"))
      confBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<CDIM> >("confBasis");
    else
      throw Lucee::Except("BGKCollUpdater::readInput: Must specify configuration space basis to use using 'confBasis'");

    if (tbl.hasNumber("mass"))
      mass = tbl.getNumber("mass");
    else
      throw Lucee::Except("BGKCollUpdater::readInput: Must specify  mass using 'mass'");

     if (tbl.hasNumber("elemCharge"))
      elemCharge = tbl.getNumber("elemCharge");
    else
      throw Lucee::Except("BGKCollUpdater::readInput: Must specify elemental charge using 'elemCharge'");

    if (tbl.hasNumber("permitivity"))
      permitivity = tbl.getNumber("permitivity");
    else
      throw Lucee::Except("BGKCollUpdater::readInput: Must specify permitivity using 'permitivity'");

    // calculate constants parts of the collision frequency and plasma
    // parameter
    collFreqConst = elemCharge*elemCharge*elemCharge*elemCharge/
      (2*M_PI*permitivity*permitivity*mass*mass);
    double temp = sqrt(permitivity*mass/(elemCharge*elemCharge));
    plasmaParamConst = temp*temp*temp;
  }
  //----------------------------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  void 
  BGKCollUpdater<CDIM, VDIM>::initialize()
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
          "BGKCollUpdater::initialize: No matching configuration space node for phase-space node ");
        lce << n;
        throw lce;
      }
    }
  }
  //----------------------------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  Lucee::UpdaterStatus 
  BGKCollUpdater<CDIM, VDIM>::update(double t)
  {
    const unsigned NDIM = CDIM+VDIM;

    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    const Lucee::Field<NDIM, double>& distf =
      this->getInp<Lucee::Field<NDIM, double> >(0);
    const Lucee::Field<CDIM, double>& zerothMoment =
      this->getInp<Lucee::Field<CDIM, double> >(1);
    const Lucee::Field<CDIM, double>& firstMoment =
      this->getInp<Lucee::Field<CDIM, double> >(2);
    const Lucee::Field<CDIM, double>& secondMoment =
      this->getInp<Lucee::Field<CDIM, double> >(3);

    Lucee::Field<NDIM, double>& rhs =
      this->getOut<Lucee::Field<NDIM, double> >(0);

    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    Lucee::Region<NDIM, double> compSpace = grid.getComputationalSpace();

    Lucee::ConstFieldPtr<double> distfPtr =
      distf.createConstPtr();
    Lucee::ConstFieldPtr<double> zerothMomentPtr =
      zerothMoment.createConstPtr();
    Lucee::ConstFieldPtr<double> firstMomentPtr =
      firstMoment.createConstPtr();
    Lucee::ConstFieldPtr<double> secondMomentPtr =
      secondMoment.createConstPtr();

    Lucee::FieldPtr<double> rhsPtr = rhs.createPtr();

    Lucee::Matrix<double> phaseNodeCoords(phaseBasis->getNumNodes(), PNC);

    rhs = 0.0; // use q to store increment initially    
    
    int idx[NDIM];
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);

    unsigned numNodesPhase = phaseBasis->getNumNodes();
    unsigned numNodesConf = confBasis->getNumNodes();
    
    std::vector<double> numDens(numNodesConf);
    std::vector<double> invNumDens(numNodesConf);

    Lucee::Matrix<double> vDrift(numNodesConf, VDIM);

    Lucee::Matrix<double> vTerm(numNodesConf, VDIM);
    Lucee::Matrix<double> invVTerm2(numNodesConf, VDIM);
    Lucee::Matrix<double> invVTerm(numNodesConf, VDIM);
    
    std::vector<double> plasmaParam(numNodesConf);
    std::vector<double> collFreq(numNodesConf);

    double n;
    double w[VDIM]; // v - u (velocity with respect to the bulk velocity)
    double invVt2[VDIM];
    double invVt[VDIM];
    
    while (seq.step()) {
      seq.fillWithIndex(idx);
      rhs.setPtr(rhsPtr, idx);

      distf.setPtr(distfPtr, idx);
      zerothMoment.setPtr(zerothMomentPtr, idx);
      firstMoment.setPtr(firstMomentPtr, idx);
      secondMoment.setPtr(secondMomentPtr, idx);
      
      phaseBasis->setIndex(idx);
      phaseBasis->getNodalCoordinates(phaseNodeCoords);
      
      for (unsigned nodeIdx = 0; nodeIdx<numNodesConf; ++nodeIdx) {
	numDens[nodeIdx] = zerothMomentPtr[nodeIdx];
	invNumDens[nodeIdx] = 1/numDens[nodeIdx];

	// add number density dependance to the collision frequency
	// and plasma parameter
	collFreq[nodeIdx] = collFreqConst*numDens[nodeIdx];
	plasmaParam[nodeIdx] = plasmaParamConst*sqrt(invNumDens[nodeIdx]);

	// gather all velocity directions
	for (unsigned dim = 0; dim<VDIM; ++dim) {
	  vDrift(nodeIdx, dim) = invNumDens[nodeIdx] *
	    firstMomentPtr[nodeIdx*VDIM + dim];
	  
	  double temp = invNumDens[nodeIdx] *
	    secondMomentPtr[nodeIdx*VDIM + dim] -
	    vDrift(nodeIdx, dim)*vDrift(nodeIdx, dim);
	  if (temp < 0)
	    return Lucee::UpdaterStatus(false, 0);
	  vTerm(nodeIdx, dim) = temp;
	  invVTerm2(nodeIdx, dim) = 1/temp;
	  invVTerm(nodeIdx, dim) = sqrt(invVTerm2(nodeIdx, dim));

	  /* THIS IS A HACK 

	     The code will only work for 1V. It is easier to adapt it
	     for 3V, but 2V is a bit tricky. It is left to be dealt
	     with in the future. */
	  collFreq[nodeIdx] *= invVTerm(nodeIdx, dim)*
	    invVTerm(nodeIdx, dim)*invVTerm(nodeIdx, dim);
	  plasmaParam[nodeIdx] *= vTerm(nodeIdx, dim)*
	    vTerm(nodeIdx, dim)*vTerm(nodeIdx, dim);
	}
	// I need to check 2 pi factor...
	collFreq[nodeIdx] *= log(2*M_PI*plasmaParam[nodeIdx]);
      }
      for (unsigned nodeIdx = 0; nodeIdx<numNodesPhase; ++nodeIdx) {
	// get number density
	n = numDens[phaseConfMap[nodeIdx]];
	
	// get volocity for Maxwellian
	for (unsigned dim = 0; dim<VDIM; ++dim)
	  w[dim] = phaseNodeCoords(nodeIdx, CDIM+dim)-
	    vDrift(phaseConfMap[nodeIdx], dim);
	
	// get thermal velocity
	for (unsigned dim = 0; dim<VDIM; ++dim) {
	  invVt2[dim] = invVTerm2(phaseConfMap[nodeIdx], dim);
	  invVt[dim] = invVTerm(phaseConfMap[nodeIdx], dim);
	}
	
	// BGK Righ-hand-side
	rhsPtr[nodeIdx] = collFreq[nodeIdx]*
	  (evaluateMaxwell(n, w, invVt, invVt2) - distfPtr[nodeIdx]);
      }
    }
    
    return Lucee::UpdaterStatus(true, 0);
  }
  //----------------------------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  void BGKCollUpdater<CDIM, VDIM>::declareTypes()
  {
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }
  //----------------------------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  double BGKCollUpdater<CDIM, VDIM>::evaluateMaxwell(double n, 
						      double w[VDIM], 
						      double invVt[VDIM],
						      double invVt2[VDIM])
  {
    double result = n;
    for (unsigned dim = 0; dim<VDIM; ++dim) {
      result *= maxwellNorm*invVt[dim];
      result *= exp(-0.5*w[dim]*w[dim]*invVt2[dim]);
    }
    return result;
  }
  //----------------------------------------------------------------------------
  // instantiations
  template class BGKCollUpdater<1, 1>;
  template class BGKCollUpdater<1, 2>;
  template class BGKCollUpdater<1, 3>;
  template class BGKCollUpdater<2, 2>;
  template class BGKCollUpdater<2, 3>;
  //template class BGKCollUpdater<3, 3>;
}

