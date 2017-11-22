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

// std includes
#include <cmath>
#include <vector>

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

    // get number of nodes in  phase space
    unsigned numNodesPhase = phaseBasis->getNumNodes();

    // get volume interpolation matrices for phase space element
    int numVolQuadPhase = phaseBasis->getNumGaussNodes();
    volWeightsPhase.resize(numVolQuadPhase);
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
    invMassMatrixPhase = massMatrixPhase.inverse();
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
    int loVel[NDIM], upVel[NDIM];

    Lucee::Field<NDIM, double>& distf =
      this->getOut<Lucee::Field<NDIM, double> >(0);
    Lucee::FieldPtr<double> distfPtrSkin = distf.createPtr();
    Lucee::FieldPtr<double> distfPtrGhost = distf.createPtr();

    // get field for node coordinates (both velocity and energetic) 
    unsigned numNodes = phaseBasis->getNumNodes();
    Lucee::Matrix<double> phaseNodeCoordsV(numNodes, PNC);  
    Lucee::Matrix<double> phaseNodeCoordsESkin(numNodes, 2);  
    Lucee::Matrix<double> phaseNodeCoordsEGhost(numNodes, 2);

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
    Lucee::Region<NDIM, int> velRgn = 
      Lucee::Region<NDIM, int>(loVel, upVel);
    // .. and sequencer
    Lucee::RowMajorSequencer<NDIM> velSeq(velRgn);

    Eigen::MatrixXd rhsMatrix(numNodes, numNodes);
    Eigen::VectorXd resultVector(numNodes);
    Eigen::VectorXd distfVector(numNodes);
    double E, cosTheta, integralResult;
    while (ghostSeq.step()) {
      ghostSeq.fillWithIndex(idxGhost);
      distf.setPtr(distfPtrGhost, idxGhost);

      ghostSeq.fillWithIndex(idxSkin);
      if (edge == LC_LOWER_EDGE)
	idxSkin[dir] = distf.getLower(dir);
      else
	idxSkin[dir] = distf.getUpper(dir);

      phaseBasis->setIndex(idxGhost);
      phaseBasis->getNodalCoordinates(phaseNodeCoordsV);
      for (int node = 0; node < numNodes; ++node) {
	E = 0;
	cosTheta = 0;
	for (int dim = 0; dim < VDIM; dim++)
	  E += phaseNodeCoordsV(node, CDIM + dim) *
	      phaseNodeCoordsV(node, CDIM + dim);
	cosTheta = sqrt((E - phaseNodeCoordsV(node, CDIM + dir) * 
			   phaseNodeCoordsV(node, CDIM + dir) )/E);
	phaseNodeCoordsEGhost(node, 0) = E;
	phaseNodeCoordsEGhost(node, 1) = cosTheta;

	distfPtrGhost[node] = 0;
      }

      while (velSeq.step()) {
	int idx[VDIM];
	velSeq.fillWithIndex(idx);
	for (unsigned d = 0; d < VDIM; ++d)
	  idxSkin[CDIM + d] = idx[d]; 
	distf.setPtr(distfPtrSkin, idxSkin);

	phaseBasis->setIndex(idxGhost);
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
	for (int i = 0; i < numNodes; ++i) {
	  for (int j = 0; j < numNodes; ++j) {
	    integralResult = 0;
	    for (int gaussIdx = 0; gaussIdx<volWeightsPhase.size(); ++gaussIdx)
	      integralResult +=
		//Rsee(phaseNodeCoordsESkin(gaussIdx, 0),
		//     phaseNodeCoordsESkin(gaussIdx, 1),
		//     phaseNodeCoordsEGhost(gaussIdx, 0),
		//     phaseNodeCoordsEGhost(gaussIdx, 1)) * 
		volWeightsPhase[gaussIdx] * 
		volQuadPhase(gaussIdx, i) * 
		volQuadPhase(gaussIdx, j);
	    rhsMatrix(i, j) = integralResult;
	  }
	}
	rhsMatrix = invMassMatrixPhase*rhsMatrix;

	for (int node = 0; node < numNodes; ++node)
	  distfVector(node) = distfPtrSkin[node];
	resultVector.noalias() = rhsMatrix*distfVector;
	
	for (int node = 0; node < numNodes; ++node)
	  distfPtrGhost[node] += resultVector(node);
      }
      velSeq.reset();
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
  //-- Copy Lucee to Eigen -------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  void
  NodalContinuumKineticSEE<CDIM, VDIM>::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); ++rowIndex)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); ++colIndex)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
 
  //------------------------------------------------------------------
  //-- Instantiations ------------------------------------------------
  template class NodalContinuumKineticSEE<1, 1>;
  template class NodalContinuumKineticSEE<1, 2>;
  template class NodalContinuumKineticSEE<1, 3>;
  template class NodalContinuumKineticSEE<2, 2>;
  template class NodalContinuumKineticSEE<2, 3>;
}

