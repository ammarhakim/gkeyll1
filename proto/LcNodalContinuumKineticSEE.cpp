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
  }

  //------------------------------------------------------------------
  //-- Update --------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  Lucee::UpdaterStatus 
  NodalContinuumKineticSEE<CDIM, VDIM>::update(double t)
  {
    const unsigned NDIM = CDIM+VDIM;
    int idxGhost[NDIM], idxSkin[NDIM];
    int loGhost[NDIM], upGhost[NDIM];
    int loVel[VDIM], upVel[VDIM];

    Lucee::Field<NDIM, double>& distf =
      this->getOut<Lucee::Field<NDIM, double> >(0);
    Lucee::FieldPtr<double> distfPtrSkin = distf.createPtr();
    Lucee::FieldPtr<double> distfPtrGhost = distf.createPtr();

    // get field for node coordinates (both velocity and energetic) 
    unsigned numNodes = phaseBasis->getNumNodes();
    Lucee::Matrix<double> phaseNodeCoords(numNodes, PNC);  
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


    double EIn, cosThetaIn, EOut, cosThetaOut;
    while (ghostSeq.step()) {
      ghostSeq.fillWithIndex(idxGhost);
      distf.setPtr(distfPtrGhost, idxGhost);

      ghostSeq.fillWithIndex(idxSkin);
      if (edge == LC_LOWER_EDGE)
	idxSkin[dir] = idxSkin[dir] + 1;
      else
	idxSkin[dir] = idxSkin[dir] - 1;

      phaseBasis->setIndex(idxGhost);
      phaseBasis->getNodalCoordinates(phaseNodeCoords);
      for (int nodeGhost = 0; nodeGhost < numSurfNodes; ++nodeGhost) {
	EOut = 0;
	for (int dim = 0; dim < VDIM; dim++)
	  EOut += phaseNodeCoords(surfNodesGhost[nodeGhost], CDIM + dim) *
	      phaseNodeCoords(surfNodesGhost[nodeGhost], CDIM + dim);
	if (EOut != 0) {
	  cosThetaOut = 
	    fabs(phaseNodeCoords(surfNodesGhost[nodeGhost], CDIM+dir)) / 
	    sqrt(EOut); 
	  EOut = 0.5*mass*EOut/elemCharge;
	} else
	  cosThetaOut = 1;

	// inner loop over skin cells
	distfPtrGhost[surfNodesGhost[nodeGhost]] = 0.0;
	while (velSeq.step()) {
	  int idx[VDIM];
	  velSeq.fillWithIndex(idx);
	  for (unsigned dim = 0; dim < VDIM; ++dim)
	    idxSkin[CDIM + dim] = idx[dim]; 
	  distf.setPtr(distfPtrSkin, idxSkin);
	  
	  phaseBasis->setIndex(idxSkin);
	  phaseBasis->getNodalCoordinates(phaseNodeCoords);

	  for (int nodeSkin = 0; nodeSkin < numSurfNodes; ++nodeSkin) {
	    EIn = 0;
	    for (int dim = 0; dim < VDIM; ++dim)
	      EIn += phaseNodeCoords(surfNodesSkin[nodeSkin], CDIM + dim) *
		phaseNodeCoords(surfNodesSkin[nodeSkin], CDIM + dim);
	    if (EIn != 0) {
	      cosThetaIn =
		fabs(phaseNodeCoords(surfNodesSkin[nodeSkin], CDIM+dir)) / 
		sqrt(EIn);
	      EIn = 0.5*mass*EIn/elemCharge;

	      distfPtrGhost[surfNodesGhost[nodeGhost]] += 
		Reflect(EIn, cosThetaIn, EOut, cosThetaOut) * 
		distfPtrSkin[surfNodesSkin[nodeSkin]] * 
		surfWeightsSkin[nodeSkin] * mass;
	    }
	  }
	}
	velSeq.reset();
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
  //-- Reflection Function -------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  double NodalContinuumKineticSEE<CDIM, VDIM>::Reflect(double EIn, 
						       double cosThetaIn, 
						       double EOut,
						       double cosThetaOut)
  {
    double delta_e0 = P_inf + (P_hat-P_inf) * 
      exp(-pow(fabs(EIn-E_hat)/W, p_e)/p_e);
    double delta_e = delta_e0 * (1 + e_1*(1 - pow(cosThetaIn, e_2)));
    double f_e = delta_e * 2 *
      exp(-(EOut-EIn)*(EOut-EIn)/(2*sigma_e*sigma_e)) / 
      (sqrt(2*M_PI)*sigma_e*erf(EIn/(sqrt(2)*sigma_e)));
    return f_e * cosThetaIn * cosThetaOut;
  }
 
  //------------------------------------------------------------------
  //-- Instantiations ------------------------------------------------
  template class NodalContinuumKineticSEE<1, 1>;
  template class NodalContinuumKineticSEE<1, 2>;
  template class NodalContinuumKineticSEE<1, 3>;
  template class NodalContinuumKineticSEE<2, 2>;
  template class NodalContinuumKineticSEE<2, 3>;
}

