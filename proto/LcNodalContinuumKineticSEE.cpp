/**
 * @file	LcNodalContinuumKineticSEE
 *
 * @brief	Updater for calculation Secondary Emitted Electron (SEE) BC
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcNodalContinuumKineticSEE.h>

// loki includes
#include <loki/Singleton.h>

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
  //template <> const char *NodalContinuumKineticSEE<3, 3>::id =
  //"NodalContinuumKineticSEE3X3V";

  //------------------------------------------------------------------
  //-- Helper Functions ----------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  bool NodalContinuumKineticSEE<CDIM,VDIM>::sameConfigCoords(unsigned n, unsigned cn, double dxMin,
							     const Lucee::Matrix<double>& phaseC,
							     const Lucee::Matrix<double>& confC)
  {
    for (unsigned dim = 0; dim < CDIM; ++dim)
      if (! (std::fabs(phaseC(n, dim)-confC(cn, dim))<1e-4*dxMin) )
        return false;
    return true;
  }
  //------------------------------------------------------------------
  //-- Reflection Function -------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  double NodalContinuumKineticSEE<CDIM, VDIM>::reflect(double eIn, double muIn, 
						       double eOut, double muOut,
						       double epsilon)
  {
    double f_e = 0;
    if (eOut - eIn < epsilon) {
      double delta_e0 = P_inf + (P_hat-P_inf) * 
	exp(-pow(fabs(eIn-E_hat)/W, p_e)/p_e);
      double delta_e = delta_e0 *
	(1 + e_1*(1 - pow(muIn, e_2)));
      f_e = delta_e * 2 *
	exp(-(eOut-eIn)*(eOut-eIn)/(2*sigma_e*sigma_e)) / 
	(sqrt(2*M_PI)*sigma_e*erf(eIn/(sqrt(2)*sigma_e)));
    }
    return f_e * fabs(muOut);
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

    if (tbl.hasNumber("polyOrder"))
      polyOrder = (unsigned) tbl.getNumber("polyOrder");
    else
      throw Lucee::Except("NodalContinuumKineticSEE::readInput: Must specify polynomial order using 'polyOrder'");

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
      if (tbl.hasNumber("dir"))
	dir = tbl.getNumber("dir");
      else
	throw Lucee::Except("NodalContinuumKineticSEE::readInput: Must specify direction using 'dir'");
    } else 
      dir = 0;

    if (tbl.hasNumber("mass")) 
      mass = tbl.getNumber("mass");
    else
      throw Lucee::Except("NodalContinuumKineticSEE::readInput: Must specify mass using 'mass'");

    if (tbl.hasNumber("elemCharge")) 
      elemCharge = tbl.getNumber("elemCharge");
    else
      throw Lucee::Except("NodalContinuumKineticSEE::readInput: Must specify elementary charge using 'elemCharge'");

    Escale = tbl.getNumber("Escale");

    // Furman model -- needs to be done properly
    E_hat = tbl.getNumber("E_hat");
    P_hat = tbl.getNumber("P_hat");
    P_inf = tbl.getNumber("P_inf");
    W = tbl.getNumber("W");
    p_e = tbl.getNumber("p_e");
    e_1 = tbl.getNumber("e_1");
    e_2 = tbl.getNumber("e_2");
    sigma_e = tbl.getNumber("sigma_e");
  }

  //------------------------------------------------------------------
  //-- Initialize ----------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  void 
  NodalContinuumKineticSEE<CDIM, VDIM>::initialize()
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
          "MaxwellDistInit::initialize: No matching configuration space node for phase-space node ");
        lce << n;
        throw lce;
      }
    }
  }

  //------------------------------------------------------------------
  //-- Update --------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  Lucee::UpdaterStatus 
  NodalContinuumKineticSEE<CDIM, VDIM>::update(double t)
  {
    const unsigned NDIM = CDIM+VDIM;
    int idxGhost[NDIM], idxSkin[NDIM];
    int idxConf[CDIM], idxVel[VDIM];
    int loGhost[NDIM], upGhost[NDIM];
    int loVel[VDIM], upVel[VDIM];

    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->template getGrid<Lucee::StructuredGridBase<NDIM> >();

    const Lucee::Field<CDIM, double>& emField =
      this->getInp<Lucee::Field<CDIM, double> >(0);
    Lucee::ConstFieldPtr<double> emFieldPtr =
      emField.createConstPtr();

    Lucee::Field<NDIM, double>& distf =
      this->getOut<Lucee::Field<NDIM, double> >(0);
    Lucee::FieldPtr<double> distfPtrSkin = distf.createPtr();
    Lucee::FieldPtr<double> distfPtrGhost = distf.createPtr();

    // get field for node coordinates (both velocity and energetic) 
    unsigned numNodes = phaseBasis->getNumNodes();
    Lucee::Matrix<double> nodeCoordsGhost(numNodes, PNC);  
    Lucee::Matrix<double> nodeCoordsSkin(numNodes, PNC);  
    // get surface nodes
    unsigned numSurfNodes = phaseBasis->getNumSurfLowerNodes(dir);
    std::vector<int> surfNodesSkin(numSurfNodes),
      surfNodesGhost(numSurfNodes);
    std::vector<double> surfWeightsSkin(numSurfNodes);
    //nodalBasis->setIndex(idx);
    if (edge == LC_LOWER_EDGE)
    {
      phaseBasis->getSurfLowerNodeNums(dir, surfNodesSkin);
      phaseBasis->getSurfLowerWeights(dir, surfWeightsSkin);
      phaseBasis->getSurfUpperNodeNums(dir, surfNodesGhost);
    }
    else
    {
      phaseBasis->getSurfUpperNodeNums(dir, surfNodesSkin);
      phaseBasis->getSurfUpperWeights(dir, surfWeightsSkin);
      phaseBasis->getSurfLowerNodeNums(dir, surfNodesGhost);
    }

    // Create a ghost layer (phase-space) region ...
    for (unsigned d = 0; d < NDIM; ++d) { 
      loGhost[d] = distf.getGlobalLowerExt(d);
      upGhost[d] = distf.getGlobalUpperExt(d);
    }
    if (edge == LC_LOWER_EDGE)
      upGhost[dir] = distf.getGlobalLower(dir);
    else
      loGhost[dir] = distf.getGlobalUpper(dir);
    Lucee::Region<NDIM, int> ghostRgn = distf.getExtRegion().intersect(
        Lucee::Region<NDIM, int>(loGhost, upGhost));
    // ... and sequencer
    Lucee::RowMajorSequencer<NDIM> ghostSeq(ghostRgn);

    // Create a velocity region ...
    for (unsigned d = 0; d < VDIM; ++d) { 
      loVel[d] = distf.getLower(CDIM + d);
      upVel[d] = distf.getUpper(CDIM + d);
    }
    Lucee::Region<VDIM, int> velRgn = 
      Lucee::Region<VDIM, int>(loVel, upVel);
    // .. and sequencer
    Lucee::RowMajorSequencer<VDIM> velSeq(velRgn);

    // analytically Ein <= Eout, however in the code we need Eout - Ein < epsilon
    double minDv = grid.getDx(CDIM);
    if (VDIM > 1)
      for (int d = 1; d < VDIM; ++d)
	minDv = std::min(minDv, grid.getDx(CDIM + d));
    double epsilon = 0.5*mass*minDv*minDv / (polyOrder+1) / elemCharge; // energy increment in eV

    double eIn, muIn, eOut, muOut, vOut, vIn, Ewall, Efactor;
    while (ghostSeq.step()) {
      ghostSeq.fillWithIndex(idxGhost);
      distf.setPtr(distfPtrGhost, idxGhost);

      ghostSeq.fillWithIndex(idxSkin);
      if (edge == LC_LOWER_EDGE)
	idxSkin[dir] = idxSkin[dir] + 1;
      else
	idxSkin[dir] = idxSkin[dir] - 1;
      for (int d = 0; d < CDIM; ++d)
	idxConf[d] = idxSkin[d];

      phaseBasis->setIndex(idxGhost);
      phaseBasis->getNodalCoordinates(nodeCoordsGhost);
      for (int n = 0; n < numNodes; ++n)
	distfPtrGhost[n] = 0.0;

      emField.setPtr(emFieldPtr, idxConf);
      Ewall = emFieldPtr[phaseConfMap[surfNodesSkin[0]]*8 + dir];
      if (edge == LC_LOWER_EDGE)
	Efactor = exp(-Ewall/Escale);
      else
	Efactor = exp(Ewall/Escale);

      for (int nodeGhost = 0; nodeGhost < numSurfNodes; ++nodeGhost) {
	eOut = 0;
	for (int dim = 0; dim < VDIM; ++dim)
	  eOut += nodeCoordsGhost(surfNodesGhost[nodeGhost], CDIM + dim) *
	      nodeCoordsGhost(surfNodesGhost[nodeGhost], CDIM + dim);
	vOut = sqrt(eOut);
	if (eOut != 0) {
	  muOut = nodeCoordsGhost(surfNodesGhost[nodeGhost], CDIM + dir) / vOut;
	  eOut = 0.5 * mass * eOut / elemCharge;
	} else
	  continue;

	if ((edge == LC_LOWER_EDGE && muOut > 0) || (edge == LC_UPPER_EDGE && muOut < 0)) {
	  // inner loop over skin cells
	  while (velSeq.step()) {
	    velSeq.fillWithIndex(idxVel);
	    for (unsigned dim = 0; dim < VDIM; ++dim)
	      idxSkin[CDIM + dim] = idxVel[dim]; 
	    distf.setPtr(distfPtrSkin, idxSkin);

	    phaseBasis->setIndex(idxSkin);
	    phaseBasis->getNodalCoordinates(nodeCoordsSkin);
	    for (int nodeSkin = 0; nodeSkin < numSurfNodes; ++nodeSkin) {
	      eIn = 0;
	      for (int dim = 0; dim < VDIM; ++dim)
		eIn += nodeCoordsSkin(surfNodesSkin[nodeSkin], CDIM + dim) *
		  nodeCoordsSkin(surfNodesSkin[nodeSkin], CDIM + dim);
	      vIn = sqrt(eIn);
	      eIn = 0.5 * mass * eIn / elemCharge;
	      if (eIn > 1) { // Artificially limit to incoming energies > 1 eV
		muIn = nodeCoordsSkin(surfNodesSkin[nodeSkin], CDIM + dir) / vIn;
		
		if ((edge == LC_UPPER_EDGE &&  muIn > 0)
		    || (edge == LC_LOWER_EDGE && muIn < 0)) {
		  double tmp = reflect(eIn, muIn, eOut, muOut, epsilon) * 
		    distfPtrSkin[surfNodesSkin[nodeSkin]] * 
		    surfWeightsSkin[nodeSkin] * mass / elemCharge * Efactor;
		  if (VDIM == 1)
		    distfPtrGhost[surfNodesGhost[nodeGhost]] += tmp * vIn;
		  else if (VDIM == 2)
		    distfPtrGhost[surfNodesGhost[nodeGhost]] += tmp;
		  else
		    distfPtrGhost[surfNodesGhost[nodeGhost]] += tmp / (2*M_PI*vIn);
		}
	      }
	    }
	  }
	  velSeq.reset();
	}
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
}

