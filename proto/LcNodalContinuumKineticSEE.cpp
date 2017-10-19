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
 
    /*const unsigned NDIM = CDIM+VDIM;

    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    
    // get handle on the local region
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[NDIM];
    seq.fillWithIndex(idx);
    phaseBasis->setIndex(idx);*/

    // get number of nodes in  phase space
    unsigned numNodesPhase = phaseBasis->getNumNodes();

    // get volume interpolation matrices for phase space element
    int numVolQuadPhase = phaseBasis->getNumGaussNodes();
    volWeightsPhase.resize(numVolQuadPhase.size());
    Lucee::Matrix<double> tempVolQuadPhase(numVolQuadPhase, numNodesPhase);
    Lucee::Matrix<double> tempVolCoordsPhase(numVolQuadPhase, (unsigned) PNC);

    phaseBasis->getGaussQuadData(tempVolQuadPhase,
				 tempVolCoordsPhase,
				 volWeightsPhase);

    volQuadPhase = Eigen::MatrixXd::Zero(numVolQuadPhase, numNodesPhase);
    copyLuceeToEigen(tempVolQuadPhase, volQuadPhase);

    // Get phase space mass matrix
    Lucee::Matrix<double> tempMassMatrixPhase(numNodesPhase, numNodesPhase);
    phaseBasis->getMassMatrix(tempMassMatrixPhase);
    Eigen::MatrixXd massMatrixPhase(numNodesPhase, numNodesPhase);
    copyLuceeToEigen(tempMassMatrixPhase, massMatrixPhase);
    invMassMatrixPhase = Eigen::MatrixXd::Zero(numNodesPhase, numNodesPhase);
    invmassMatrixPhase = massMatrixPhase.inverse()
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

    Lucee::Field<NDIM, double>& distf =
      this->getOut<Lucee::Field<NDIM, double> >(0);
    Lucee::FieldPtr<double> distfPtrSkin = distf.createPtr();
    Lucee::FieldPtr<double> distfPtrGhost = distf.createPtr();

    // get field for node coordinates (velocity and energetic) 
    unsigned numNodes = phaseBasis->getNumNodes();
    Lucee::Matrix<double> phaseNodeCoordsV(numNodes, PNC);  
    Lucee::Matrix<double> phaseNodeCoordsESkin(numNodes, 2);  
    Lucee::Matrix<double> phaseNodeCoordsEGhost(numNodes, 2);

    // Creating regions
    Lucee::Region<NDIM, int> localRgn = distf.getRegion();
    Lucee::Region<NDIM, int> externalRgn = distf.getExtRegion();

    int lowerConfSkin[CDIM], upperConfSkin[CDIM];
    int lowerConfGhost[CDIM], upperConfGhost[CDIM];
    for (int i = 0; i < CDIM; i++) {
      lowerConfEdge[i] = localRgn.getLower(i);
      upperConfEdge[i] = localRgn.getUpper(i);
      lowerConfGhost[i] = localRgn.getLower(i);
      upperConfGhost[i] = localRgn.getUpper(i);
    }
    if (edge == LC_UPPER_EDGE) {
      lowerConfEdge[dir] = localRng.getUpper(dir);
      lowerConfGhost[dir] = externalRng.getUpper(dir);
      upperConfGhost[dir] = externalRng.getUpper(dir);
    } else {
      upperConfEdge[dir] = localRng.getLower(dir);
      lowerConfGhost[dir] = externalRng.getLower(dir);
      upperConfGhost[dir] = externalRng.getLower(dir);
    }
    Lucee::Region<CDIM, int> confSkinRgn(lowerConfSkin, upperConfSkin);
    Lucee::Region<CDIM, int> confGhostRgn(lowerConfGhost, upperConfGhost);

    int lowerVel[CDIM], upperVel[CDIM];
    for (int i = 0; i < VDIM; i++) {
      lowerVel[i] = localRgn.getLower(CDIM + i);
      upperVel[i] = localRgn.getUpper(CDIM + i);
    }
    Lucee::Region<VDIM, int> velRgn(lowerVel, upperVel);
    
    // Creating sequencers
    int idxSkin[NDIM], idxGhost[NDIM], idxVel[VDIM];
    Lucee::RowMajorSequencer<CDIM> seqConfSkin(confSkinRgn);
    Lucee::RowMajorSequencer<CDIM> seqConfGhost(confGhostRgn);
    Lucee::RowMajorSequencer<VDIM> seqVelSkin(velRgn);
    Lucee::RowMajorSequencer<VDIM> seqVelGhost(velRgn);

    // Edge nodes
    unsigned numSurfNodes = phaseBasis->getNumSurfLowerNodes(dir);
    std::vector<int> surfNodesSkin(numSurfNodes), 
      surfNodesGhost(numSurfNodes);
    if (edge == LC_LOWER_EDGE) {
      phaseBasis->getSurfUpperNodeNums(dir, surfNodesGhost);
      phaseBasis->getSurfLowerNodeNums(dir, surfNodesSkin);
    } else {
      phaseBasis->getSurfLowerNodeNums(dir, surfNodesGhost);
      phaseBasis->getSurfUpperNodeNums(dir, surfNodesSkin);
    }

    Eigen::MatrixXd rhsMatrix(numNodes, numNodes);
    Eigen::VectorXd resultVector(numNodes), distfVector(numNodes);
    double E, cosTheta, integralResult;
    while (seqConfSkin.step()) {
      seqConfSkin.fillWithIndex(idxSkin);
      seqConfGhost.step();
      seqConfGhost.fillWithIndex(idxGhost);

      while (seqVelGhost.step()) {
	seqVelGhost.fillWithIndex(*idxGhost + CDIM);
	distf.setPtr(distfPtrGhost, idxGhost);

	phaseBasis->setIndex(idxGhost);
	phaseBasis->getNodalCoordinates(phaseNodeCoordsV);

	// convert coordinates to energies and angles
	for (int node = 0; node < numNodes; ++node) {
	  E = 0;
	  cosTheta = 0;
	  for (int dim = 0; dim < VDIM; dim++)
	    E += phaseNodeCoordsV(node, CDIM + dim) *
	      phaseNodeCoordsV(node, CDIM + dim);
	  cosTheta = sqrt((E - phaseNodeCoordsV(node, CDIM+dir) * 
			   phaseNodeCoordsV(node, CDIM+dir))/E);
	  phaseNodeCoordsEGhost(node, 0) = E;
	  phaseNodeCoordsEGhost(node, 1) = cosTheta;

	  distfPtrGhost[node] = 0;
	}

	while (seqVelSkin.step()) {
	  seqVelSkin.fillWithIndex(*idxSkin + CDIM);
	  distf.setPtr(distfPtrSkin, idxSkin);
	  
	  // convert coordinates to energies and angles
	  phaseBasis->setIndex(idxSkin);
	  phaseBasis->getNodalCoordinates(phaseNodeCoordsV);
	  for (int node = 0; node < numNodes; ++node) {
	    E = 0;
	    cosTheta = 0;
	    for (int dim = 0; dim < VDIM; dim++)
	      E += phaseNodeCoordsV(node, CDIM + dim) *
		phaseNodeCoordsV(node, CDIM + dim);
	    cosTheta = sqrt((E - phaseNodeCoordsV(node, CDIM+dir) * 
			     phaseNodeCoordsV(node, CDIM+dir))/E);
	    phaseNodeCoordsESkin(node, 0) = E;
	    phaseNodeCoordsESkin(node, 1) = cosTheta;
	  }

	  // integrate
	  for (int i = 0; i < numNodes, ++i) {
	    for (int j = 0; j < numNodes, ++j) {
	      integralResult = 0;
	      for (int gaussIdx = 0; gaussIdx < volWeightsPhase.size(); ++gaussIdx)
		integralResult +=
		  getSEE(phaseNodeCoordsESkin(gaussIdx, 0),
			 phaseNodeCoordsESkin(gaussIdx, 1),
			 phaseNodeCoordsEGhost(gaussIdx, 0),
			 phaseNodeCoordsEGhost(gaussIdx, 1)) * 
		  volWeightsPhase[gaussIdx] * 
		  volQuadPhase(gaussIdx, i) * 
		  volQuadPhase(gaussIdx, j);
	      rhsMatrix(i, j) = integralResult;
	    }
	  }
	  rhsMatrix = invMassMatrixPhase*rhsMatrix;

	  for (int node = 0; node < numNodes, ++node)
	    distfVector(node) = distfPtrSkin[node];
	  resultVector.noalias() = rhsMatrix*distfVector;

	  for (int node = 0; node < numSurfNodes, ++node)
	    distfPtrGhost[surfNodesGhost[node]] = resultVector(surfNodesSkin[node])
	  
	}
	seqVelIn.reset();
      }
      seqVelOut.reset();
    }

    while (seqOut.step()) {
      seq.fillWithIndex(idx);
      q.setPtr(qPtr, idx);
      zerothMoment.setPtr(zerothMomentPtr, idx);
      firstMoment.setPtr(firstMomentPtr, idx);
      secondMoment.setPtr(secondMomentPtr, idx);
      

      
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

