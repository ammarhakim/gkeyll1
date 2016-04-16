/**
 * @file	LcBGKCollUpdater.cpp
 *
 * @brief	Updater to get a right-hand-side for the Boltzmann eqation. For more details see "Greene; Improved Bhatnagar-Gross-Krook model of electron-ion collisions; The Physics of fluids; 1973".
 */

// lucee includes
#include <LcBGKCollUpdater.h>

namespace Lucee
{
// set ids for module system
  template <> const char *BGKCollUpdater<1, 1>::id =
    "BGKCollUpdater1X1V";
  template <> const char *BGKCollUpdater<1, 2>::id =
    "BGKCollUpdater1X2V";
  template <> const char *BGKCollUpdater<1, 3>::id =
    "BGKCollUpdater1X3V";
  template <> const char *BGKCollUpdater<2, 2>::id =
    "BGKCollUpdater2X2V";
  template <> const char *BGKCollUpdater<2, 3>::id =
    "BGKCollUpdater2X3V";
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

    // Collision frequency within the species itself
    if (tbl.hasNumber("nuSelf"))
      nuSelf = tbl.getNumber("nuSelf");
    else
      throw Lucee::Except("BGKCollUpdater::readInput: Must specify a collision frequency within the species itself to use using 'nuSelf'");

    // Collision frequency with the other species
    if (tbl.hasNumber("nuCross"))
      nuCross = tbl.getNumber("nuCross");
    else
      throw Lucee::Except("BGKCollUpdater::readInput: Must specify a collision frequency with the other species to use using 'nuCross'");

    // Mass of the species
    if (tbl.hasNumber("massSelf"))
      massSelf = tbl.getNumber("massSelf");
    else
      throw Lucee::Except("BGKCollUpdater::readInput: Must specify a mass of the species to use using 'massSelf'");

    // Mass of the other species
    if (tbl.hasNumber("massOther"))
      massSelf = tbl.getNumber("massOther");
    else
      throw Lucee::Except("BGKCollUpdater::readInput: Must specify a mass of the other species to use using 'massOther'");

    // Arbitrary parameter beta
    if (tbl.hasNumber("beta"))
      beta = tbl.getNumber("beta");
    else
      throw Lucee::Except("BGKCollUpdater::readInput: Must specify an arbitrary parameter beta to use using 'beta'. See 'Greene; Improved Bhatnagar-Gross-Krook model of electron-ion collisions; The Physics of fluids; 1973' for details.");
    
    invMassSelf = 1/massSelf;
    invMassOther = 1/massOther;
    invMassSum = 1/(massSelf+massOther);
    precalculation = (1-beta*beta)*massSelf*massOther/6+(1+beta*beta)*massSelf*(massOther-massSelf)/12;
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
    const double sixth = 1.0/6.0;
    const double twelfth = 1.0/12.0;

    const unsigned NDIM = CDIM+VDIM;

    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    // Inputs
    const Lucee::Field<NDIM, double>& distFn =
      this->getInp<Lucee::Field<NDIM, double> >(0);
    const Lucee::Field<CDIM, double>& zerothMomentSelf =
      this->getInp<Lucee::Field<CDIM, double> >(1);
    const Lucee::Field<CDIM, double>& firstMomentSelf =
      this->getInp<Lucee::Field<CDIM, double> >(2);
    const Lucee::Field<CDIM, double>& secondMomentSelf =
      this->getInp<Lucee::Field<CDIM, double> >(3);
    const Lucee::Field<CDIM, double>& zerothMomentOther =
      this->getInp<Lucee::Field<CDIM, double> >(4);
    const Lucee::Field<CDIM, double>& firstMomentOther =
      this->getInp<Lucee::Field<CDIM, double> >(5);
    const Lucee::Field<CDIM, double>& secondMomentOther =
      this->getInp<Lucee::Field<CDIM, double> >(6);
 
    // Output
    // For compatibility with mathematical definition of DG, I am
    // retaining the names "q" for the distribution function.
    Lucee::Field<NDIM, double>& q =
      this->getOut<Lucee::Field<NDIM, double> >(0);

    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    Lucee::Region<NDIM, double> compSpace = grid.getComputationalSpace();

    // Creating updaters to inputs and outputs
    Lucee::FieldPtr<double> qPtr = q.createPtr();
    Lucee::ConstFieldPtr<double> distFnPtr = 
      distFn.createConstPtr();
    Lucee::ConstFieldPtr<double> zerothMomentSelfPtr =
      zerothMomentSelf.createConstPtr();
    Lucee::ConstFieldPtr<double> firstMomentSelfPtr =
      firstMomentSelf.createConstPtr();
    Lucee::ConstFieldPtr<double> secondMomentSelfPtr =
      secondMomentSelf.createConstPtr();
    Lucee::ConstFieldPtr<double> zerothMomentOtherPtr =
      zerothMomentOther.createConstPtr();
    Lucee::ConstFieldPtr<double> firstMomentOtherPtr =
      firstMomentOther.createConstPtr();
    Lucee::ConstFieldPtr<double> secondMomentOtherPtr =
      secondMomentOther.createConstPtr();

    Lucee::Matrix<double> phaseNodeCoords(phaseBasis->getNumNodes(), PNC);

    q = 0.0; // use q to store increment initially    
    
    int idx[NDIM];
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);

    unsigned numNodesPhase = phaseBasis->getNumNodes();
    unsigned numNodesConf = confBasis->getNumNodes();
    
    // Maxwellian parameters to be calculated from moments
    std::vector<double> densSelf(numNodesConf);
    std::vector<double> densOther(numNodesConf);
    std::vector<double> invDensSelf(numNodesConf);
    std::vector<double> invDensOther(numNodesConf);
    std::vector<double> vTerm2Self(numNodesConf);
    std::vector<double> vTerm2Other(numNodesConf);
    Lucee::Matrix<double> vDriftSelf(numNodesConf, VDIM);
    Lucee::Matrix<double> vDriftOther(numNodesConf, VDIM);
    double vSelf[VDIM];
    double vCross[VDIM];
    double densSelfMapped;
    double vTerm2SelfMapped;

    // Variables needed to calculate the "cross Maxwellian"
    std::vector<double> TSelf(numNodesConf);
    std::vector<double> TOther(numNodesConf);
    std::vector<double> vTerm2Cross(numNodesConf);
    double vTerm2CrossMapped;
    std::vector<double> vDriftDiff2(numNodesConf);
    double vDriftDiff;
    Lucee::Matrix<double> vDriftCross(numNodesConf, VDIM);

    // get an "area" of velocity space so the updater can return a
    // uniform distribution with correct density in case the higher
    // moments were not calculated correctly and vTerm2 < 0
    double vSpaceArea = 1;
    double vSpaceAreaInv;
    for (int dim = 0; dim<VDIM; ++dim)
    {
      vSpaceArea *= compSpace.getShape(CDIM+dim);
    }
    vSpaceAreaInv = 1/vSpaceArea;

    while (seq.step())
      {
	// Set up the pointers
	seq.fillWithIndex(idx);
	q.setPtr(qPtr, idx);
	distFn.setPtr(distFnPtr, idx);
	zerothMomentSelf.setPtr(zerothMomentSelfPtr, idx);
	firstMomentSelf.setPtr(firstMomentSelfPtr, idx);
	secondMomentSelf.setPtr(secondMomentSelfPtr, idx);
	zerothMomentOther.setPtr(zerothMomentOtherPtr, idx);
	firstMomentOther.setPtr(firstMomentOtherPtr, idx);
	secondMomentOther.setPtr(secondMomentOtherPtr, idx);

	phaseBasis->setIndex(idx);
	phaseBasis->getNodalCoordinates(phaseNodeCoords);
      
	// Calculated drift and thermal velocity from moments
	for (unsigned nodeIdx = 0; nodeIdx<numNodesConf; ++nodeIdx)
	{
	  densSelf[nodeIdx] = zerothMomentSelfPtr[nodeIdx];
	  invDensSelf[nodeIdx] = 1/zerothMomentSelfPtr[nodeIdx];
	  vTerm2Self[nodeIdx] = invDensSelf[nodeIdx]*secondMomentSelfPtr[nodeIdx];
	  densOther[nodeIdx] = zerothMomentOtherPtr[nodeIdx];
	  invDensOther[nodeIdx] = 1/zerothMomentOtherPtr[nodeIdx];
	  vTerm2Other[nodeIdx] = invDensOther[nodeIdx]*secondMomentOtherPtr[nodeIdx];
	  for (unsigned dim = 0; dim<VDIM; ++dim)
	  {
	    vDriftSelf(nodeIdx, dim) = invDensSelf[nodeIdx]*
	      firstMomentSelfPtr[nodeIdx*VDIM+dim];
	    vTerm2Self[nodeIdx] = vTerm2Self[nodeIdx]-
	      vDriftSelf(nodeIdx, dim)*vDriftSelf(nodeIdx, dim);
	    vDriftOther(nodeIdx, dim) = invDensOther[nodeIdx]*
	      firstMomentOtherPtr[nodeIdx*VDIM+dim];
	    vTerm2Other[nodeIdx] = vTerm2Other[nodeIdx]-
	      vDriftOther(nodeIdx, dim)*vDriftOther(nodeIdx, dim);

	    TSelf[nodeIdx] = vTerm2Self[nodeIdx]*massSelf;
	    TOther[nodeIdx] = vTerm2Other[nodeIdx]*massOther;
	  }  	
	}

	// Calculate drift and thermal velocity of the "cross Maxwellian"
	for (unsigned nodeIdx = 0; nodeIdx<numNodesConf; ++nodeIdx)
	{
	  vDriftDiff2[nodeIdx] = 0;
	  for (unsigned dim = 0; dim<VDIM; ++dim)
	  {
	    vDriftDiff = vDriftSelf(nodeIdx, dim)-vDriftOther(nodeIdx, dim);
	    vDriftCross(nodeIdx, dim) = 
	      0.5*(vDriftSelf(nodeIdx, dim)+vDriftOther(nodeIdx, dim)) -
	      0.5*beta*vDriftDiff;
	    vDriftDiff2[nodeIdx] = vDriftDiff2[nodeIdx]+vDriftDiff*vDriftDiff;
	  }  
	  vTerm2Cross[nodeIdx] = (massOther*TSelf[nodeIdx]+massSelf*TOther[nodeIdx] - 
				  beta*massSelf*(TSelf[nodeIdx]-TOther[nodeIdx]) + 
				  vDriftDiff2[nodeIdx]*precalculation)*invMassSelf*invMassSum;
	  /** precalculation = (1-beta*beta)*massSelf*massOther/6+
	      (1+beta*beta)*massSelf*(massother-massSelf)/12; 
	      This value is calculated during the initialization of the updater*/
	}

	for (unsigned nodeIdx = 0; nodeIdx<numNodesPhase; ++nodeIdx) 
	{
	  for (unsigned dim = 0; dim<VDIM; ++dim)
	  {
	    vSelf[dim] = phaseNodeCoords(nodeIdx, CDIM+dim)-
	      vDriftSelf(phaseConfMap[nodeIdx], dim);
	    vCross[dim] = phaseNodeCoords(nodeIdx, CDIM+dim)-
	      vDriftCross(phaseConfMap[nodeIdx], dim);
	  }
	  densSelfMapped = densSelf[phaseConfMap[nodeIdx]];
	  vTerm2SelfMapped = vTerm2Self[phaseConfMap[nodeIdx]];
	  vTerm2CrossMapped = vTerm2Cross[phaseConfMap[nodeIdx]];
	  //if (vTerm2In >= 0) 
	  //  qPtr[nodeIdx] = evaluateMaxwell(v, densIn, vTerm2In);
	  //else
	  //  qPtr[nodeIdx] = densIn*vSpaceAreaInv;
	  qPtr[nodeIdx] = 
	    nuSelf*(evaluateMaxwell(vSelf, densSelfMapped, vTerm2SelfMapped)-distFnPtr[nodeIdx])+
	    nuCross*(evaluateMaxwell(vCross, densSelfMapped, vTerm2CrossMapped)-distFnPtr[nodeIdx]);
	}
      }

    return Lucee::UpdaterStatus();
  }
  //----------------------------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  void BGKCollUpdater<CDIM, VDIM>::declareTypes()
  {
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }
  //----------------------------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  double BGKCollUpdater<CDIM, VDIM>::evaluateMaxwell(double v[VDIM], double n, 
    double vTerm2)
  {
    double result = 0;
    double factor = 1/sqrt(2*M_PI*vTerm2);
    for (unsigned dim = 0; dim<VDIM; ++dim)
      result += v[dim]*v[dim];
    result = n*exp(-0.5*result/vTerm2);
    for (unsigned dim = 0; dim<VDIM; ++dim)
      result *= factor;
    return result;
  }
  //----------------------------------------------------------------------------
  // instantiations
  template class BGKCollUpdater<1, 1>;
  template class BGKCollUpdater<1, 2>;
  template class BGKCollUpdater<1, 3>;
  template class BGKCollUpdater<2, 2>;
  template class BGKCollUpdater<2, 3>;
  //template class MaxwellDistInit<3, 3>;
}

