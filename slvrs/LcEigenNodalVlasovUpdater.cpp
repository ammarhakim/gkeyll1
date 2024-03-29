/**
 * @file	LcEigenNodalVlasovUpdater.cpp
 *
 * @brief	Updater to solve Vlasov equation with nodal DG scheme using Eigen library.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAlignedRectCoordSys.h>
#include <LcField.h>
#include <LcLinAlgebra.h>
#include <LcMathLib.h>
#include <LcEigenNodalVlasovUpdater.h>
#include <LcStructuredGridBase.h>

// std includes
#include <vector>
#include <cmath>

namespace Lucee
{
// set id for module system
  template <> const char *EigenNodalVlasovUpdater<1,1>::id = "EigenNodalVlasov1X1V";
  template <> const char *EigenNodalVlasovUpdater<1,2>::id = "EigenNodalVlasov1X2V";
  template <> const char *EigenNodalVlasovUpdater<1,3>::id = "EigenNodalVlasov1X3V";
  template <> const char *EigenNodalVlasovUpdater<2,2>::id = "EigenNodalVlasov2X2V";
  template <> const char *EigenNodalVlasovUpdater<2,3>::id = "EigenNodalVlasov2X3V";
  //template <> const char *EigenNodalVlasovUpdater<3,3>::id = "EigenNodalVlasov3X3V";

// makes indexing a little more sane
  static const unsigned IX = 0;
  static const unsigned IY = 1;
  static const unsigned IZ = 2;

  static const unsigned IEX = 0;
  static const unsigned IEY = 1;
  static const unsigned IEZ = 2;
  static const unsigned IBX = 3;
  static const unsigned IBY = 4;
  static const unsigned IBZ = 5;

  // helper to index EM fields at nodes
  template <unsigned CDIM, unsigned VDIM>  
  unsigned
  EigenNodalVlasovUpdater<CDIM,VDIM>::emidx(unsigned n, unsigned i)
  {
    return n*8+i; // 8 as last two are correction potentials
  }

  // helpers that returns 0 if not enough velocity space dimensions
  template <unsigned CDIM, unsigned VDIM>
  double
  EigenNodalVlasovUpdater<CDIM,VDIM>::getSafeVx(int n, const Lucee::Matrix<double>& pc)
  { return pc(n, CDIM+IX); }

  template <unsigned CDIM, unsigned VDIM>  
  double
  EigenNodalVlasovUpdater<CDIM,VDIM>::getSafeVy(int n, const Lucee::Matrix<double>& pc)
  { return VDIM>1 ? pc(n, CDIM+IY) : 0; }

  template <unsigned CDIM, unsigned VDIM>
  double
  EigenNodalVlasovUpdater<CDIM,VDIM>::getSafeVz(int n, const Lucee::Matrix<double>& pc)
  { return VDIM>2 ? pc(n, CDIM+IZ) : 0; }

  template <unsigned CDIM, unsigned VDIM>
  bool
  EigenNodalVlasovUpdater<CDIM,VDIM>::sameConfigCoords(unsigned n, unsigned cn, double dxMin,
    const Eigen::MatrixXd& phaseC, const Eigen::MatrixXd& confC)
  {
    for (unsigned d=0; d<CDIM; ++d)
      if (! (std::fabs(phaseC(n,d)-confC(cn,d))<1e-4*dxMin) )
        return false;
    return true;
  }

  template <unsigned CDIM, unsigned VDIM>
  EigenNodalVlasovUpdater<CDIM,VDIM>::EigenNodalVlasovUpdater()
    : UpdaterIfc()
  {
  }

  template <unsigned CDIM, unsigned VDIM>  
  void 
  EigenNodalVlasovUpdater<CDIM,VDIM>::readInput(Lucee::LuaTable& tbl)
  {
    const unsigned NDIM = CDIM+VDIM;
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("phaseBasis"))
      phaseBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("phaseBasis");
    else
      throw Lucee::Except("EigenNodalVlasovUpdater::readInput: Must specify phase-space basis using 'phaseBasis'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<CDIM> >("confBasis"))
      confBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<CDIM> >("confBasis");
    else
      throw Lucee::Except("EigenNodalVlasovUpdater::readInput: Must specify configuration-space basis using 'confBasis'");

    if (tbl.hasNumber("polyOrder"))
      polyOrder = tbl.getNumber("polyOrder");
    else
      throw Lucee::Except("EigenNodalVlasovUpdater::readInput: Must specify polynomial order with 'polyOrder' for anisotropic quadrature");

    skipVelocitySweep = false;
    if (tbl.hasBool("skipVelocitySweep"))
      skipVelocitySweep = tbl.getBool("skipVelocitySweep");

    applyZeroFluxBc = true;
    if (tbl.hasBool("applyZeroFluxBc"))
      applyZeroFluxBc = tbl.getBool("applyZeroFluxBc");

    cfl = tbl.getNumber("cfl");
    cflm = 1.1*cfl; // use slightly large max CFL to avoid thrashing around

    onlyIncrement = false;
// when onlyIncrement flag is set contribution is not added to the
// input field, i.e. only increment is computed
    if (tbl.hasBool("onlyIncrement"))
      onlyIncrement = tbl.getBool("onlyIncrement");

    charge = tbl.getNumber("charge");
    mass = tbl.getNumber("mass");
  }

  template <unsigned CDIM, unsigned VDIM>
  void 
  EigenNodalVlasovUpdater<CDIM,VDIM>::initialize()
  {
    const unsigned NDIM = CDIM+VDIM;
    Lucee::UpdaterIfc::initialize();

// get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
// local region to update
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[NDIM];
    seq.fillWithIndex(idx);
    phaseBasis->setIndex(idx);
    confBasis->setIndex(idx); // only first CDIM elements are used
    
    unsigned nlocal = phaseBasis->getNumNodes();

// get node numbers on each lower and upper edges since only those nodes needed for surface updates
    for (unsigned dir=0; dir<NDIM; ++dir)
    {
      lowerNodeNums[dir].nums.resize(phaseBasis->getNumSurfLowerNodes(dir));
      phaseBasis->getSurfLowerNodeNums(dir, lowerNodeNums[dir].nums);

      upperNodeNums[dir].nums.resize(phaseBasis->getNumSurfUpperNodes(dir));
      phaseBasis->getSurfUpperNodeNums(dir, upperNodeNums[dir].nums);
    }

// compute mapping of phase-space nodes to configuration space
// nodes. The assumption here is that the node layout in phase-space
// and configuration space are such that each node in phase-space has
// exactly one node co-located with it in configuration space. No
// "orphan" phase-space node are allowed, and an exception is thrown
// if that occurs.
    phaseConfMap.resize(nlocal);
    Lucee::Matrix<double> tempPhaseNodeCoords(phaseBasis->getNumNodes(), PNC);
    Lucee::Matrix<double> tempConfNodeCoords(confBasis->getNumNodes(), CNC);

    phaseBasis->getNodalCoordinates(tempPhaseNodeCoords);
    confBasis->getNodalCoordinates(tempConfNodeCoords);

    Eigen::MatrixXd phaseNodeCoords(phaseBasis->getNumNodes(), (unsigned) PNC);
    copyLuceeToEigen(tempPhaseNodeCoords, phaseNodeCoords);

    Eigen::MatrixXd confNodeCoords(confBasis->getNumNodes(), (unsigned) CNC);
    copyLuceeToEigen(tempConfNodeCoords, confNodeCoords);

    double dxMin = grid.getDx(0);
    for (unsigned d=1; d<CDIM; ++d)
      dxMin = std::min(dxMin, grid.getDx(d));

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
          "EigenNodalVlasovUpdater::readInput: No matching configuration space node for phase-space node ");
        lce << n;
        throw lce;
      }
    }

    // Store mass matrix inverse
    Lucee::Matrix<double> tempMass(nlocal, nlocal);
    phaseBasis->getMassMatrix(tempMass);
    Eigen::MatrixXd massMatrix(nlocal, nlocal);
    copyLuceeToEigen(tempMass, massMatrix);
    massMatrixInv = massMatrix.inverse();

    // Store grad stiffness matrix in each direction
    std::vector<Eigen::MatrixXd> gradStiffnessMatrix(NDIM);
    for (int dir=0; dir<NDIM; ++dir)
    {
      Lucee::Matrix<double> tempMatrix(nlocal, nlocal);
      phaseBasis->getGradStiffnessMatrix(dir, tempMatrix);
      gradStiffnessMatrix[dir] = Eigen::MatrixXd(nlocal, nlocal);

      copyLuceeToEigen(tempMatrix, gradStiffnessMatrix[dir]);
    }

    //Number of quadrature points required for streaming term and velocity space in forcing term
    int numGaussPoints1d = ceil((2*polyOrder+1)/2.0);
    std::vector<double> gaussPoints1d(numGaussPoints1d);
    std::vector<double> gaussWeights1d(numGaussPoints1d);
    legendre_set(numGaussPoints1d, &gaussPoints1d[0], &gaussWeights1d[0]);

    //Number of quadrature points required for configuration space in forcing term
    int numGaussPoints1dField = ceil((3*polyOrder+1)/2.0);
    std::vector<double> gaussPoints1dField(numGaussPoints1dField);
    std::vector<double> gaussWeights1dField(numGaussPoints1dField);
    legendre_set(numGaussPoints1dField, &gaussPoints1dField[0], &gaussWeights1dField[0]);

    //Temporary array to evaluate basis functions at quadrature points to
    std::vector<double> evalBasisTemp(nlocal);
    //Temporary array to store ordinates of quadrature points
    double xcTemp[NDIM];

    double weightScale = 1.0;
    for (int dimIndex=0; dimIndex<NDIM; ++dimIndex)
      weightScale *= 0.5*grid.getDx(dimIndex);

    // Figure out how many gaussian points there are
    nvolQuadStream = 1;
    nvolQuadForce = 1;

    for (int dimIndex=0; dimIndex<NDIM; ++dimIndex)
    {
      nvolQuadStream *= gaussPoints1d.size();
      if (dimIndex < CDIM)
      {
        nvolQuadForce *= gaussPoints1dField.size();
      }
      else
      {
        nvolQuadForce *= gaussPoints1d.size();
      }
    }

    volQuadStream.reset(nvolQuadStream, nlocal, PNC);
    volQuadForce.reset(nvolQuadForce, nlocal, PNC);

    Eigen::MatrixXd gaussNodeListStream = Eigen::MatrixXd::Zero(nvolQuadStream, NDIM+1);
    Eigen::MatrixXd gaussNodeListForce = Eigen::MatrixXd::Zero(nvolQuadForce, NDIM+1);

    // Evaluate all basis functions at all gaussian volume nodes
    // First compute all volume gaussian quadrature locations
    int volShapeStream[NDIM];
    int volShapeForce[NDIM];
    int volIdxStream[NDIM];
    int volIdxForce[NDIM];
    
    for (int dimIndex=0; dimIndex<NDIM; ++dimIndex)
    {
      volShapeStream[dimIndex] = gaussPoints1d.size();
      if (dimIndex < CDIM)
      {
        volShapeForce[dimIndex] = gaussPoints1dField.size();
      }
      else
      {
        volShapeForce[dimIndex] = gaussPoints1d.size();
      }
    }

    Lucee::Region<NDIM, int> volRegionStream(volShapeStream);
    Lucee::ColMajorSequencer<NDIM> volSeqStream = ColMajorSequencer<NDIM>(volRegionStream);
    Lucee::ColMajorIndexer<NDIM> volIdxrStream = ColMajorIndexer<NDIM>(volRegionStream);

    // Find all quadrature locations on a NDIM volume for streaming term
    while(volSeqStream.step())
    {
      volSeqStream.fillWithIndex(volIdxStream);
      int nodeNumberStream = volIdxrStream.getIndex(volIdxStream);

      gaussNodeListStream(nodeNumberStream, NDIM) = 1.0;
      for (int dimIndex=0; dimIndex<NDIM; ++dimIndex)
      {
        xcTemp[dimIndex] = gaussPoints1d[volIdxStream[dimIndex]];
        gaussNodeListStream(nodeNumberStream, dimIndex) = gaussPoints1d[volIdxStream[dimIndex]];
        gaussNodeListStream(nodeNumberStream, NDIM)    *= gaussWeights1d[volIdxStream[dimIndex]];
      }
      volQuadStream.weights(nodeNumberStream) = weightScale*gaussNodeListStream(nodeNumberStream, NDIM);
      phaseBasis->evalBasis(xcTemp, evalBasisTemp);
      for(int basisIndex=0; basisIndex<nlocal; ++basisIndex)
      {
        volQuadStream.interpMat(nodeNumberStream, basisIndex) = evalBasisTemp[basisIndex];
      }
    }

    Lucee::Region<NDIM, int> volRegionForce(volShapeForce);
    Lucee::ColMajorSequencer<NDIM> volSeqForce = ColMajorSequencer<NDIM>(volRegionForce);
    Lucee::ColMajorIndexer<NDIM> volIdxrForce = ColMajorIndexer<NDIM>(volRegionForce);

    // Find all quadrature locations on a NDIM volume for force term
    while(volSeqForce.step())
    {
      volSeqForce.fillWithIndex(volIdxForce);
      int nodeNumberForce = volIdxrForce.getIndex(volIdxForce);

      gaussNodeListForce(nodeNumberForce, NDIM) = 1.0;
      for (int dimIndex=0; dimIndex<NDIM; ++dimIndex)
      {
        if (dimIndex < CDIM)
        {
          xcTemp[dimIndex] = gaussPoints1dField[volIdxForce[dimIndex]];
          gaussNodeListForce(nodeNumberForce, dimIndex) = gaussPoints1dField[volIdxForce[dimIndex]];
          gaussNodeListForce(nodeNumberForce, NDIM)    *= gaussWeights1dField[volIdxForce[dimIndex]];
        }
        else
        {
          xcTemp[dimIndex] = gaussPoints1d[volIdxForce[dimIndex]];
          gaussNodeListForce(nodeNumberForce, dimIndex) = gaussPoints1d[volIdxForce[dimIndex]];
          gaussNodeListForce(nodeNumberForce, NDIM)    *= gaussWeights1d[volIdxForce[dimIndex]];
        }
      }
      volQuadForce.weights(nodeNumberForce) = weightScale*gaussNodeListForce(nodeNumberForce, NDIM);
      phaseBasis->evalBasis(xcTemp, evalBasisTemp);
      for(int basisIndex=0; basisIndex<nlocal; ++basisIndex)
      {
        volQuadForce.interpMat(nodeNumberForce, basisIndex) = evalBasisTemp[basisIndex];
      }
    }

    std::vector<Eigen::MatrixXd> derivMatrices(VDIM);
    Eigen::MatrixXd derivMatrix;

    // Compute gradients of basis functions evaluated at volume quadrature points
    for (int dir=0; dir<VDIM; ++dir)
    {
      int updateDir = dir+CDIM;
      // Each row is a quadrature point; each column is a basis function with derivative applied
      derivMatrix.noalias() = volQuadForce.interpMat*massMatrixInv*gradStiffnessMatrix[updateDir].transpose();
      derivMatrices[dir] = derivMatrix;
    }

    /* Since the streaming term just involves a coordinate, we can generate a quadrature free update by 
       employing a local velocity coordinate, exactly like how we evaluate the moments in the moment 
       calculation routine */

    volStreamFluxMatrices = std::vector<Eigen::MatrixXd>(CDIM); 
    volGradStiffnessMatrices = std::vector<Eigen::MatrixXd>(CDIM);
    for (int h=0; h<CDIM; ++h)
    {
      volStreamFluxMatrices[h] = Eigen::MatrixXd::Zero(nlocal, nlocal);
      volGradStiffnessMatrices[h] = Eigen::MatrixXd::Zero(nlocal, nlocal);
      int updateDir = h+CDIM;
      for (int i=0; i<nlocal; ++i)
      {
        for (int j=0; j<nlocal; ++j)
        {
          double integralResultCoordinate = 0.0;
          for (int gaussIndex=0; gaussIndex<nvolQuadStream; ++gaussIndex)
          {
            double baseIntegral
              = volQuadStream.weights(gaussIndex)*volQuadStream.interpMat(gaussIndex,i)*volQuadStream.interpMat(gaussIndex,j);
            // get coordinate of quadrature point in direction updateDir
            double coord2Val = gaussNodeListStream(gaussIndex, updateDir)*grid.getDx(updateDir)/2.0;
            integralResultCoordinate += coord2Val*baseIntegral;
          }
          volStreamFluxMatrices[h](i,j) = integralResultCoordinate;
        }
      }
      // multiply matrices by inverse of mass matrix
      volStreamFluxMatrices[h] = massMatrixInv*volStreamFluxMatrices[h];
      volGradStiffnessMatrices[h] = massMatrixInv*gradStiffnessMatrix[h];
    }


    // Figure out how many surface quadrature gaussian points there are
    nSurfQuadStream = 1;
    nSurfQuadForce = 1;

    for (int dimIndex=0; dimIndex<NDIM-1; ++dimIndex)
    {
      nSurfQuadStream *= gaussPoints1d.size();
      if (dimIndex < CDIM)
      {
        nSurfQuadForce *= gaussPoints1dField.size();
      }
      else
      {
        nSurfQuadForce *= gaussPoints1d.size();
      }
    }

    /* Due to the symmetry of the system, there are still only two node lists, one for the streaming term
       and one for the forcing term. This is because the anisotropic quadrature is only used for the 
       forcing term, which involves gradients in velocity space. In other words, the velocity space surfaces
       require the quadrature information of all of the corresponding position in configuration space */

    Eigen::MatrixXd gaussNodeListSurfStream = Eigen::MatrixXd::Zero(nSurfQuadStream, NDIM);
    Eigen::MatrixXd gaussNodeListSurfForce = Eigen::MatrixXd::Zero(nSurfQuadForce, NDIM);

    int surfShapeStream[NDIM-1];
    int surfIdxStream[NDIM-1];
    int surfShapeForce[NDIM-1];
    int surfIdxForce[NDIM-1];

    for (int dimIndex=0; dimIndex<NDIM-1; ++dimIndex)
    {
      surfShapeStream[dimIndex] = gaussPoints1d.size();
      if (dimIndex < CDIM)
      {
        surfShapeForce[dimIndex] = gaussPoints1dField.size();
      }
      else
      {
        surfShapeForce[dimIndex] = gaussPoints1d.size();
      }
    }

    Lucee::Region<NDIM-1, int> surfRegionStream(surfShapeStream);
    Lucee::ColMajorSequencer<NDIM-1> surfSeqStream = ColMajorSequencer<NDIM-1>(surfRegionStream);
    Lucee::ColMajorIndexer<NDIM-1> surfIdxrStream = ColMajorIndexer<NDIM-1>(surfRegionStream);

    // Find all quadrature locations on a NDIM-1 surface for streaming term
    while(surfSeqStream.step())
    {
      surfSeqStream.fillWithIndex(surfIdxStream);
      int nodeNumberStream = surfIdxrStream.getIndex(surfIdxStream);

      gaussNodeListSurfStream(nodeNumberStream, NDIM-1) = 1.0;
      for (int dimIndex=0; dimIndex<NDIM-1; ++dimIndex)
      {
        gaussNodeListSurfStream(nodeNumberStream, dimIndex) = gaussPoints1d[surfIdxStream[dimIndex]];
        gaussNodeListSurfStream(nodeNumberStream, NDIM-1)    *= gaussWeights1d[surfIdxStream[dimIndex]];
      }
    }

    Lucee::Region<NDIM-1, int> surfRegionForce(surfShapeForce);
    Lucee::ColMajorSequencer<NDIM-1> surfSeqForce = ColMajorSequencer<NDIM-1>(surfRegionForce);
    Lucee::ColMajorIndexer<NDIM-1> surfIdxrForce = ColMajorIndexer<NDIM-1>(surfRegionForce);


    // Find all quadrature locations on a NDIM-1 surface for force term
    while(surfSeqForce.step())
    {
      surfSeqForce.fillWithIndex(surfIdxForce);
      int nodeNumberForce = surfIdxrForce.getIndex(surfIdxForce);

      gaussNodeListSurfForce(nodeNumberForce, NDIM-1) = 1.0;
      for (int dimIndex=0; dimIndex<NDIM-1; ++dimIndex)
      {
        if (dimIndex < CDIM)
        {
          gaussNodeListSurfForce(nodeNumberForce, dimIndex) = gaussPoints1dField[surfIdxForce[dimIndex]];
          gaussNodeListSurfForce(nodeNumberForce, NDIM-1)    *= gaussWeights1dField[surfIdxForce[dimIndex]];
        }
        else
        {
          gaussNodeListSurfForce(nodeNumberForce, dimIndex) = gaussPoints1d[surfIdxForce[dimIndex]];
          gaussNodeListSurfForce(nodeNumberForce, NDIM-1)    *= gaussWeights1d[surfIdxForce[dimIndex]];
        }
      }
    }

    // Get data for surface quadrature
    for (int dir=0; dir<VDIM; ++dir)
    {

      double weightScale = 1.0;
      int updateDir = dir + CDIM;
      for (int dimIndex=0; dimIndex<NDIM; ++dimIndex)
      {
        if (dimIndex != updateDir)
          weightScale *= 0.5*grid.getDx(dimIndex);
      }

      unsigned nface = lowerNodeNums[updateDir].nums.size();

      // Reset surface quadrature structures
      surfFullLowerQuad[dir].reset(nSurfQuadForce, nlocal, PNC);
      surfFullUpperQuad[dir].reset(nSurfQuadForce, nlocal, PNC);
      surfLowerQuad[dir].reset(nSurfQuadForce, nface, PNC);
      surfUpperQuad[dir].reset(nSurfQuadForce, nface, PNC);

      // Evaluate all basis functions at all upper and lower surfaces
      for (int nodeIndex = 0; nodeIndex < gaussNodeListSurfForce.rows(); nodeIndex++)
      {
        int rollingIndex = 0;
        
        for (int coordIndex = 0; coordIndex < NDIM; coordIndex++)
        {
          if (coordIndex == updateDir)
            xcTemp[updateDir] = 1;
          else
          {
              xcTemp[coordIndex] = gaussNodeListSurfForce(nodeIndex, rollingIndex);
              rollingIndex++;
          }
        }
        //Upper and lower surface quadrature weights are the same, so compute and then copy in
        surfFullUpperQuad[dir].weights(nodeIndex) = weightScale*gaussNodeListSurfForce(nodeIndex, NDIM-1);
        surfFullLowerQuad[dir].weights(nodeIndex) = surfFullUpperQuad[dir].weights(nodeIndex);

        //Evaluate basis functions at upper surface quadrature nodes
        phaseBasis->evalBasis(xcTemp, evalBasisTemp);
        for(int basisIndex=0; basisIndex<nlocal; ++basisIndex)
        {
          surfFullUpperQuad[dir].interpMat(nodeIndex, basisIndex) = evalBasisTemp[basisIndex];
        }
        for(int basisIndex=0; basisIndex<nface; ++basisIndex)
        {
          unsigned un = upperNodeNums[updateDir].nums[basisIndex];
          surfUpperQuad[dir].interpMat(nodeIndex, basisIndex) = evalBasisTemp[un];
        }

        // Replace surface index for computation for lower surface matrices
        xcTemp[updateDir] = -1;
        phaseBasis->evalBasis(xcTemp, evalBasisTemp);
        for(int basisIndex=0; basisIndex<nlocal; ++basisIndex)
        {
          surfFullLowerQuad[dir].interpMat(nodeIndex, basisIndex) = evalBasisTemp[basisIndex];
        }
        for(int basisIndex=0; basisIndex<nface; ++basisIndex)
        {
          unsigned ln = lowerNodeNums[updateDir].nums[basisIndex];
          surfLowerQuad[dir].interpMat(nodeIndex, basisIndex) = evalBasisTemp[ln];
        }
      }
    }

    /* Similar to the volume integral, since the streaming term just involves a coordinate, 
       we can generate a quadrature free update by employing a local velocity coordinate;
       however, here we must take care to keep track of the local coordinate on both the 
       upper and lower surface */

    //These matrices are used to transform to the lower coordinate
    upperSurfStreamFluxMatrices = std::vector<Eigen::MatrixXd>(CDIM); 
    lowerSurfStreamFluxMatrices = std::vector<Eigen::MatrixXd>(CDIM);
    surfPermutationUpperMatrices = std::vector<Eigen::MatrixXd>(CDIM); 
    surfPermutationLowerMatrices = std::vector<Eigen::MatrixXd>(CDIM);
    //These matrices are identical to the upper and lower mass matrices 
    //and then multiplied by the inverse of the mass matrix
    upperLiftMatrices = std::vector<Eigen::MatrixXd>(CDIM);
    lowerLiftMatrices = std::vector<Eigen::MatrixXd>(CDIM);

    //Temporary arrays to hold lift matrices before copying them in to Eigen arrays
    Lucee::Matrix<double> lowerLift;
    Lucee::Matrix<double> upperLift;

    for (int h=0; h<CDIM; ++h)
    {
      unsigned nface = lowerNodeNums[h].nums.size();
      // compute lift matrices
      upperLift = Lucee::Matrix<double>(nlocal, nface);

      phaseBasis->getUpperFaceMassMatrix(h, upperLift);

      lowerLift = Lucee::Matrix<double>(nlocal, nface);

      phaseBasis->getLowerFaceMassMatrix(h, lowerLift);

      upperSurfStreamFluxMatrices[h] = Eigen::MatrixXd::Zero(nface, nlocal);
      lowerSurfStreamFluxMatrices[h] = Eigen::MatrixXd::Zero(nface, nlocal);

      surfPermutationUpperMatrices[h] = Eigen::MatrixXd::Zero(nface, nlocal);
      surfPermutationLowerMatrices[h] = Eigen::MatrixXd::Zero(nface, nlocal);

      upperLiftMatrices[h] = Eigen::MatrixXd::Zero(nlocal, nface);
      lowerLiftMatrices[h] = Eigen::MatrixXd::Zero(nlocal, nface);

      for (int j=0; j<nface; ++j)
      {
        unsigned un = upperNodeNums[h].nums[j];
        unsigned ln = lowerNodeNums[h].nums[j];

        surfPermutationUpperMatrices[h](j, un) = 1.0;
        surfPermutationLowerMatrices[h](j, ln) = 1.0; 
      }

      copyLuceeToEigen(upperLift, upperLiftMatrices[h]);
      copyLuceeToEigen(lowerLift, lowerLiftMatrices[h]);

      // multiply matrices by inverse of mass matrix
      upperSurfStreamFluxMatrices[h] = surfPermutationUpperMatrices[h]*volStreamFluxMatrices[h];
      lowerSurfStreamFluxMatrices[h] = surfPermutationLowerMatrices[h]*volStreamFluxMatrices[h];

      upperLiftMatrices[h] = massMatrixInv*upperLiftMatrices[h];
      lowerLiftMatrices[h] = massMatrixInv*lowerLiftMatrices[h];
    }

    bigStoredUpperSurfMatrices.resize(VDIM);
    bigStoredLowerSurfMatrices.resize(VDIM);
    bigStoredVolMatrices.resize(VDIM);

    // Store three matrices at each cell
    for (int dir=0; dir<VDIM; ++dir)
    {
      bigStoredUpperSurfMatrices[dir] = massMatrixInv*surfFullUpperQuad[dir].interpMat.transpose();
      bigStoredLowerSurfMatrices[dir] = massMatrixInv*surfFullLowerQuad[dir].interpMat.transpose();
      bigStoredVolMatrices[dir] = massMatrixInv*derivMatrices[dir].transpose();
      for (int i=0; i<nlocal; ++i)
      {
        bigStoredUpperSurfMatrices[dir].row(i) = bigStoredUpperSurfMatrices[dir].row(i).cwiseProduct(surfFullUpperQuad[dir].weights.transpose());
        bigStoredLowerSurfMatrices[dir].row(i) = bigStoredLowerSurfMatrices[dir].row(i).cwiseProduct(surfFullLowerQuad[dir].weights.transpose());
        bigStoredVolMatrices[dir].row(i) = bigStoredVolMatrices[dir].row(i).cwiseProduct(volQuadForce.weights.transpose());
      }
    }

    if (applyZeroFluxBc)
    {
// initialize directions in which zero-flux BCs are applied
      for (unsigned d=0; d<CDIM; ++d)
        lowerZeroFluxOffset[d] = upperZeroFluxOffset[d] = 0; // NO at configuration-space edges
      for (unsigned d=CDIM; d<NDIM; ++d)
        lowerZeroFluxOffset[d] = upperZeroFluxOffset[d] = 1; // YES at velocity-space edges

// ensure that zero-flux BCs are applied only if local rank owns the
// skin cell    
      Lucee::Region<NDIM, int> globalRgn = grid.getGlobalRegion();
      for (unsigned d=CDIM; d<NDIM; ++d)
      {
        if (localRgn.getLower(d) != globalRgn.getLower(d))
          lowerZeroFluxOffset[d] = 0; // not owned by us, so ignore
        if (localRgn.getUpper(d) != globalRgn.getUpper(d))
          upperZeroFluxOffset[d] = 0; // not owned by us, so ignore
      }
    }
    else
    {
// requested NOT to apply zero-flux BCs
      for (unsigned d=0; d<NDIM; ++d)
        lowerZeroFluxOffset[d] = upperZeroFluxOffset[d] = 0;
    }

    tm1 = tm2 = 0.0;
  }

  template <unsigned CDIM, unsigned VDIM>
  Lucee::UpdaterStatus 
  EigenNodalVlasovUpdater<CDIM,VDIM>::update(double t)
  {
    const unsigned NDIM = CDIM+VDIM;    
    const Lucee::StructuredGridBase<NDIM>& grid
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

// for compatibility with NodalDisContHyperUpdater I am retaining the
// names "q" and "qNew" for the distribution function. Hence, q is the
// distribution function at time t and qNew the distribution function
// at t+dt (Ammar Hakim)
    const Lucee::Field<NDIM, double>& q = this->getInp<Lucee::Field<NDIM, double> >(0);
    const Lucee::Field<CDIM, double>& EM = this->getInp<Lucee::Field<CDIM, double> >(1);
    Lucee::Field<NDIM, double>& qNew = this->getOut<Lucee::Field<NDIM, double> >(0);

    unsigned nlocal = phaseBasis->getNumNodes();
    int nSurfQuad = phaseBasis->getNumSurfGaussNodes();

    double dt = t-this->getCurrTime();
    double cfla = 0.0; // maximum CFL number used
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    Lucee::ConstFieldPtr<double> qPtr = q.createConstPtr();
    Lucee::ConstFieldPtr<double> qPtrl = q.createConstPtr();
    Lucee::FieldPtr<double> qNewPtr = qNew.createPtr();
    Lucee::FieldPtr<double> qNewPtrl = qNew.createPtr();
    Lucee::ConstFieldPtr<double> emPtr = EM.createConstPtr();
    Lucee::ConstFieldPtr<double> emPtrl = EM.createConstPtr();

    Eigen::VectorXd flux(nlocal);
    Eigen::VectorXd fVec(nlocal);
    Eigen::VectorXd fStreamVec(nlocal);

    Eigen::VectorXd alphaForce(nvolQuadForce);
    Eigen::VectorXd fAtQuadForce(nvolQuadForce);

    Eigen::VectorXd resultVectorDir = Eigen::VectorXd::Zero(nlocal);
    Eigen::VectorXd resultVectorDirLeft = Eigen::VectorXd::Zero(nlocal);

    Eigen::VectorXd rightDataAtQuadForce(nSurfQuadForce);
    Eigen::VectorXd leftDataAtQuadForce(nSurfQuadForce);
    Eigen::VectorXd alphaRightForce(nSurfQuadForce);
    Eigen::VectorXd alphaLeftForce(nSurfQuadForce);
    Eigen::VectorXd maxFluxForce(nSurfQuadForce);
    Eigen::VectorXd numericalFluxAtQuadForce(nSurfQuadForce);

    Lucee::Matrix<double> phaseNodeCoords(phaseBasis->getNumNodes(), PNC);

    qNew = 0.0; // use qNew to store increment initially    
    
    int idx[NDIM];
    double xc[PNC];
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);

    clock_t t1, t2;
    t1 = clock();
// loop to compute contribution from volume integrals
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      grid.setIndex(idx);
      grid.getCentroid(xc);

      q.setPtr(qPtr, idx);
      EM.setPtr(emPtr, idx); // only CDIM indices will be used
      qNew.setPtr(qNewPtr, idx);

      phaseBasis->setIndex(idx);
      phaseBasis->getNodalCoordinates(phaseNodeCoords);

      // Get a vector of f at quad points
      for (int i=0; i<nlocal; ++i)
        fVec(i) = qPtr[i];
      fAtQuadForce.noalias() = volQuadForce.interpMat*fVec;

      for (int dir=0;  dir<NDIM; ++dir)
      {
        if (dir<CDIM)
        {
          int updateDir = dir+CDIM;
          fStreamVec.noalias() = volStreamFluxMatrices[dir]*fVec;
          fStreamVec = fStreamVec + xc[updateDir]*fVec;
          resultVectorDir.noalias() = volGradStiffnessMatrices[dir]*fStreamVec;
        }
        else
        {
          int updateVelDir = dir-CDIM;
          calcFlux(dir, phaseNodeCoords, emPtr, flux);
          alphaForce.noalias() = volQuadForce.interpMat*flux;
          resultVectorDir.noalias() = bigStoredVolMatrices[updateVelDir]*(fAtQuadForce.cwiseProduct(alphaForce));
        }
        for (int i=0; i<nlocal; ++i)
          qNewPtr[i] += resultVectorDir(i);
      }

    }
    t2 = clock();
    tm1 += (double) (t2-t1)/CLOCKS_PER_SEC;

    t1 = clock();
// contributions from surface integrals
    for (unsigned dir=0; dir<NDIM; ++dir)
    {
// create sequencer to loop over *each* 1D slice in 'dir' direction
      Lucee::RowMajorSequencer<NDIM> seq(localRgn.deflate(dir));

// lower and upper bounds of 1D slice. (We need to make sure that flux
// is computed for one edge outside domain interior, accounting for
// the fact that we may not want to compute fluxes from the outermost
// edges [zero-flux BCs])
      int sliceLower = localRgn.getLower(dir)+lowerZeroFluxOffset[dir];
      int sliceUpper = localRgn.getUpper(dir)+1-upperZeroFluxOffset[dir];

      int idx[NDIM], idxl[NDIM];
      double xc[PNC];
// loop over each 1D slice
      while (seq.step())
      {
        seq.fillWithIndex(idx);
        seq.fillWithIndex(idxl);

        for (int sliceIndex=sliceLower; sliceIndex<sliceUpper; ++sliceIndex)
        { // loop over each edge
          idx[dir] = sliceIndex; // cell right of edge
          idxl[dir] = sliceIndex-1; // cell left of edge

// set to cell right of edge (which means we need to use *left* face
// node numbers below to get coordinates)
          phaseBasis->setIndex(idx);
          phaseBasis->getNodalCoordinates(phaseNodeCoords);

          grid.setIndex(idxl);
          double dxL = grid.getDx(dir);
          grid.setIndex(idx);
          double dxR = grid.getDx(dir);
          grid.getCentroid(xc);

          double dtdx = 2*dt/(dxL+dxR);

          q.setPtr(qPtr, idx);
          q.setPtr(qPtrl, idxl);
          EM.setPtr(emPtr, idx);
          EM.setPtr(emPtrl, idxl);

// compute numerical fluxes on each face
          unsigned nface = lowerNodeNums[dir].nums.size();

          Eigen::VectorXd dataOnFace(nface);
          Eigen::VectorXd fluxOnFace(nface);

          if (dir < CDIM)
          {          
            int updateDir = dir+CDIM;

            for (int i=0; i<nlocal; ++i)
              fVec(i) = qPtr[i];
            dataOnFace.noalias() = (lowerSurfStreamFluxMatrices[dir] + xc[updateDir]*surfPermutationLowerMatrices[dir])*fVec;

            for (int i=0; i<nlocal; ++i)
              fVec(i) = qPtrl[i];
            fluxOnFace.noalias() = (upperSurfStreamFluxMatrices[dir] + xc[updateDir]*surfPermutationUpperMatrices[dir])*fVec;

            for (int i=0; i<nface; ++i)
            {
              unsigned ln = lowerNodeNums[dir].nums[i];
              cfla = Lucee::max3(cfla, dtdx*phaseNodeCoords(ln, updateDir), -dtdx*phaseNodeCoords(ln, updateDir));
              // this choice of velocity avoids problems when a node lies exactly at 0.0
              double vel = phaseNodeCoords(ln, updateDir) + 1e-4*(xc[updateDir] - phaseNodeCoords(ln, updateDir));
              if (vel > 0)
                dataOnFace(i) = fluxOnFace(i);
            }
            resultVectorDirLeft.noalias() = upperLiftMatrices[dir]*dataOnFace;
            resultVectorDir.noalias() = lowerLiftMatrices[dir]*dataOnFace;
          }
          else
          {
            int updateVelDir = dir-CDIM;
            calcFlux(dir, phaseNodeCoords, emPtr, flux);
            // Copy data to Eigen vectors
            for (int i=0; i<nface; ++i)
            {
              unsigned ln = lowerNodeNums[dir].nums[i];
              dataOnFace(i) = qPtr[ln];
              fluxOnFace(i) = flux[ln];
            }
            rightDataAtQuadForce.noalias() = surfLowerQuad[updateVelDir].interpMat*dataOnFace;
            alphaRightForce.noalias() = surfLowerQuad[updateVelDir].interpMat*fluxOnFace;

            calcFlux(dir, phaseNodeCoords, emPtrl, flux);
            // Copy data to Eigen vectors
            for (int i=0; i<nface; ++i)
            {
              unsigned un = upperNodeNums[dir].nums[i];
              dataOnFace(i) = qPtrl[un];
              fluxOnFace(i) = flux[un];
            }
            leftDataAtQuadForce.noalias() = surfUpperQuad[updateVelDir].interpMat*dataOnFace;
            alphaLeftForce.noalias() = surfUpperQuad[updateVelDir].interpMat*fluxOnFace;

            maxFluxForce = (alphaLeftForce.cwiseAbs()).cwiseMax(alphaRightForce.cwiseAbs());

            // Compute numerical flux
            numericalFluxAtQuadForce.noalias() = 0.5*(alphaLeftForce.cwiseProduct(leftDataAtQuadForce) + alphaRightForce.cwiseProduct(rightDataAtQuadForce)) 
              - 0.5*maxFluxForce.cwiseProduct(rightDataAtQuadForce - leftDataAtQuadForce);

            cfla = Lucee::max3(cfla, dtdx*maxFluxForce.maxCoeff(), -dtdx*maxFluxForce.maxCoeff());

            resultVectorDirLeft.noalias() = bigStoredUpperSurfMatrices[updateVelDir]*numericalFluxAtQuadForce;
            resultVectorDir.noalias() = bigStoredLowerSurfMatrices[updateVelDir]*numericalFluxAtQuadForce;
          }
          qNew.setPtr(qNewPtr, idx);
          qNew.setPtr(qNewPtrl, idxl);

          for (int i=0; i<nlocal; ++i)
          {
            qNewPtrl[i] -= resultVectorDirLeft(i);
            qNewPtr[i] += resultVectorDir(i);
          }
        }
      }
      // Check to see if we need to retake time step
      if (cfla > cflm)
        return Lucee::UpdaterStatus(false, dt*cfl/cfla);
    }
    t2 = clock();
    tm2 += double (t2-t1)/CLOCKS_PER_SEC;

// NOTE: If only calculation of increments are requested, the final
// Euler update is not performed. This means that the multiplication
// of the DG RHS with dt is not done, something to keep in mind if
// using the increment in time-dependent update.
    if (onlyIncrement == false)
    {
      seq = Lucee::RowMajorSequencer<NDIM>(localRgn);
// final sweep, update solution with forward Euler step
      while (seq.step())
      {
        seq.fillWithIndex(idx);
        qNew.setPtr(qNewPtr, idx);
        q.setPtr(qPtr, idx);
        for (unsigned k=0; k<nlocal; ++k)
          qNewPtr[k] = qPtr[k] + dt*qNewPtr[k];
      }
    }

    return Lucee::UpdaterStatus(true, dt*cfl/cfla);
  }

  template <unsigned CDIM, unsigned VDIM>  
  void
  EigenNodalVlasovUpdater<CDIM,VDIM>::declareTypes()
  {
    const unsigned NDIM = CDIM+VDIM;    
// distribution function
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
// E and B field in a single field
    this->appendInpVarType(typeid(Lucee::Field<CDIM, double>));
// returns one output: updated distribution function
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned CDIM, unsigned VDIM>
  void
  EigenNodalVlasovUpdater<CDIM,VDIM>::calcFlux(unsigned dir, const Lucee::Matrix<double>& pc,
    const Lucee::ConstFieldPtr<double>& EM, Eigen::VectorXd& flux)
  {
    unsigned nlocal = flux.size();
    double qbym = charge/mass;
    if (dir<CDIM)
    { // configuration space flux
      for (unsigned n=0; n<nlocal; ++n)
        flux[n] = pc(n,CDIM+dir);
    }
    else if (dir==(CDIM+0))
    { // VX flux
      for (unsigned n=0; n<nlocal; ++n)
      {
        double vy = getSafeVy(n,pc), vz = getSafeVz(n,pc);
        double Ex = EM[emidx(phaseConfMap[n],IEX)];
        double Bz = EM[emidx(phaseConfMap[n],IBZ)];
        double By = EM[emidx(phaseConfMap[n],IBY)];
        flux[n] = qbym*(Ex + vy*Bz-vz*By);
      }
    }
    else if (dir==(CDIM+1))
    { // VY flux
      for (unsigned n=0; n<nlocal; ++n)
      {
        double vx = getSafeVx(n,pc), vz = getSafeVz(n,pc);
        double Ey = EM[emidx(phaseConfMap[n],IEY)];
        double Bx = EM[emidx(phaseConfMap[n],IBX)];
        double Bz = EM[emidx(phaseConfMap[n],IBZ)];
        flux[n] = qbym*(Ey + vz*Bx-vx*Bz);
      }
    }
    else if (dir==(CDIM+2))
    { // VZ flux
      for (unsigned n=0; n<nlocal; ++n)
      {
        double vx = getSafeVx(n,pc), vy = getSafeVy(n,pc);
        double Ez = EM[emidx(phaseConfMap[n],IEZ)];
        double By = EM[emidx(phaseConfMap[n],IBY)];
        double Bx = EM[emidx(phaseConfMap[n],IBX)];
        flux[n] = qbym*(Ez + vx*By-vy*Bx);
      }
    }
  }

  template <unsigned CDIM, unsigned VDIM>
  void
  EigenNodalVlasovUpdater<CDIM,VDIM>::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex=0; rowIndex<destinationMatrix.rows(); ++rowIndex)
      for (int colIndex=0; colIndex<destinationMatrix.cols(); ++colIndex)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }

  template <unsigned CDIM, unsigned VDIM>
  void
  EigenNodalVlasovUpdater<CDIM,VDIM>::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
    lfm.appendFunc("timers", luaGetTimers);
  }  

  template <unsigned CDIM, unsigned VDIM>
  int
  EigenNodalVlasovUpdater<CDIM,VDIM>::luaGetTimers(lua_State *L)
  {
    EigenNodalVlasovUpdater<CDIM,VDIM> *s
      = Lucee::PointerHolder<EigenNodalVlasovUpdater<CDIM,VDIM> >::getObj(L);
    lua_pushnumber(L, s->tm1);
    lua_pushnumber(L, s->tm2);    
    return 2;
  }    
  
// instantiations
  template class EigenNodalVlasovUpdater<1,1>;
  template class EigenNodalVlasovUpdater<1,2>;
  template class EigenNodalVlasovUpdater<1,3>;
  template class EigenNodalVlasovUpdater<2,2>;
  template class EigenNodalVlasovUpdater<2,3>;
  //template class EigenNodalVlasovUpdater<3,3>;
}
