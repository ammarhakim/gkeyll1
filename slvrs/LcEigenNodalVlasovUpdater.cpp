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

    // Store mass matrix inverse
    Lucee::Matrix<double> tempMass(nlocal, nlocal);
    phaseBasis->getMassMatrix(tempMass);
    Eigen::MatrixXd massMatrix(nlocal, nlocal);
    copyLuceeToEigen(tempMass, massMatrix);
    massMatrixInv = massMatrix.inverse();

    // Store grad stiffness matrix in each direction
    std::vector<Eigen::MatrixXd> gradStiffnessMatrix(NDIM);
    for (int dir = 0; dir < NDIM; dir++)
    {
      Lucee::Matrix<double> tempMatrix(nlocal, nlocal);
      phaseBasis->getGradStiffnessMatrix(dir, tempMatrix);
      gradStiffnessMatrix[dir] = Eigen::MatrixXd(nlocal, nlocal);

      copyLuceeToEigen(tempMatrix, gradStiffnessMatrix[dir]);
    }

    // get number of surface quadrature points
    int nSurfQuad = phaseBasis->getNumSurfGaussNodes();
    // get data needed for Gaussian quadrature
    int nVolQuad = phaseBasis->getNumGaussNodes();
    std::vector<double> volWeights(nVolQuad);
    Lucee::Matrix<double> tempVolQuad(nVolQuad, nlocal);
    Lucee::Matrix<double> tempVolCoords(nVolQuad, PNC);
    volQuad.reset(nVolQuad, nlocal, PNC);

    phaseBasis->getGaussQuadData(tempVolQuad, tempVolCoords, volWeights);
    for (int volIndex = 0; volIndex < nVolQuad; volIndex++)
      volQuad.weights(volIndex) = volWeights[volIndex];
    
    copyLuceeToEigen(tempVolQuad, volQuad.interpMat);

    std::vector<Eigen::MatrixXd> derivMatrices(NDIM);

    // Compute gradients of basis functions evaluated at volume quadrature points
    for (int dir = 0; dir < NDIM; dir++)
    {
      // Each row is a quadrature point; each column is a basis function with derivative applied
      Eigen::MatrixXd derivMatrix = volQuad.interpMat*massMatrixInv*gradStiffnessMatrix[dir].transpose();

      derivMatrices[dir] = derivMatrix;
    }

    // Get data for surface quadrature
    for (int dir = 0; dir < NDIM; dir++)
    {
      // temporary variables
      std::vector<double> tempSurfWeights(nSurfQuad);
      Lucee::Matrix<double> tempSurfQuad(nSurfQuad, nlocal);
      Lucee::Matrix<double> tempSurfCoords(nSurfQuad, PNC);

      // Reset surface quadrature structures
      surfLowerQuad[dir].reset(nSurfQuad, nlocal, PNC);
      surfUpperQuad[dir].reset(nSurfQuad, nlocal, PNC);
      
      // lower surface data
      phaseBasis->getSurfLowerGaussQuadData(dir, tempSurfQuad,
        tempSurfCoords, tempSurfWeights);
      // copy data to appropriate structures
      for (int quadIndex = 0; quadIndex < nSurfQuad; quadIndex++)
        surfLowerQuad[dir].weights(quadIndex) = tempSurfWeights[quadIndex];
      copyLuceeToEigen(tempSurfQuad, surfLowerQuad[dir].interpMat);

      // upper surface data
      phaseBasis->getSurfUpperGaussQuadData(dir, tempSurfQuad,
        tempSurfCoords, tempSurfWeights);
      // copy data to appropriate structures
      for (int quadIndex = 0; quadIndex < nSurfQuad; quadIndex++)
        surfUpperQuad[dir].weights(quadIndex) = tempSurfWeights[quadIndex];
      copyLuceeToEigen(tempSurfQuad, surfUpperQuad[dir].interpMat);
    }

    // CONSIDER: instead of allocating NDIM size vectors, only store those needed in updateDirs
    bigStoredUpperSurfMatrices.resize(NDIM);
    bigStoredLowerSurfMatrices.resize(NDIM);
    bigStoredVolMatrices.resize(NDIM);

    // Store three matrices at each cell
    for (int dir = 0; dir < NDIM; dir++)
    {
      bigStoredUpperSurfMatrices[dir] = massMatrixInv*surfUpperQuad[dir].interpMat.transpose();
      bigStoredLowerSurfMatrices[dir] = massMatrixInv*surfLowerQuad[dir].interpMat.transpose();
      bigStoredVolMatrices[dir] = massMatrixInv*derivMatrices[dir].transpose();
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
    int nVolQuad = phaseBasis->getNumGaussNodes();
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

    Eigen::MatrixXd alpha(NDIM, nVolQuad);
    Eigen::VectorXd fAtQuad(nVolQuad);

    std::vector<Eigen::VectorXd> resultVectorDir = std::vector<Eigen::VectorXd>(NDIM);
    for (int i=0; i<NDIM; ++i)
      resultVectorDir[i] = Eigen::VectorXd::Zero(nlocal);
    std::vector<Eigen::VectorXd> resultVectorDirLeft = std::vector<Eigen::VectorXd>(NDIM);
    for (int i=0; i<NDIM; ++i)
      resultVectorDirLeft[i] = Eigen::VectorXd::Zero(nlocal);    

    Eigen::VectorXd rightData(nlocal);
    Eigen::VectorXd leftData(nlocal);
    Eigen::VectorXd rightDataAtQuad(nSurfQuad);
    Eigen::VectorXd leftDataAtQuad(nSurfQuad);
    Eigen::VectorXd alphaRight(nSurfQuad);
    Eigen::VectorXd alphaLeft(nSurfQuad);
    Eigen::VectorXd numericalFluxAtQuad(nSurfQuad);

    double localQ, localQl, localF;

    Lucee::Matrix<double> phaseNodeCoords(phaseBasis->getNumNodes(), PNC);

    qNew = 0.0; // use qNew to store increment initially    
    
    int idx[NDIM];
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);

    clock_t t1, t2;
    t1 = clock();
// loop to compute contribution from volume integrals
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      q.setPtr(qPtr, idx);
      EM.setPtr(emPtr, idx); // only CDIM indices will be used
      qNew.setPtr(qNewPtr, idx);

      phaseBasis->setIndex(idx);
      phaseBasis->getNodalCoordinates(phaseNodeCoords);

      // Get alpha from appropriate function
      for (unsigned dir=0; dir<NDIM; ++dir)
      {
        calcFlux(dir, phaseNodeCoords, emPtr, flux);
        alpha.row(dir).noalias() = volQuad.interpMat*flux;
      }

      // Computing maximum cfla
      grid.setIndex(idx);
      Eigen::VectorXd dtdqVec(NDIM);
      for (int i = 0; i < NDIM; i++)
        dtdqVec(i) = dt/grid.getDx(i);

      Eigen::VectorXd cflaVec(NDIM);
      // Loop over each alpha vector evaluated at a quadrature point
      for (int i = 0; i < nVolQuad; i++)
      {
        double maxCflaInCol = alpha.col(i).cwiseAbs().cwiseProduct(dtdqVec).maxCoeff();
        if (maxCflaInCol > cfla)
          cfla = maxCflaInCol;
      }

      // Get a vector of f at quad points
      for (int i = 0; i < nlocal; i++)
        fVec(i) = qPtr[i];
      fAtQuad.noalias() = volQuad.interpMat*fVec;

      for (int dir = 0;  dir < NDIM; dir++)
      {
        resultVectorDir[dir].noalias() = bigStoredVolMatrices[dir]*(volQuad.weights.cwiseProduct(fAtQuad.cwiseProduct(alpha.row(dir).transpose())));
      }

      for (int i = 0; i < nlocal; i++)
      {
        for (int dir = 0;  dir < NDIM; dir++)
          qNewPtr[i] += resultVectorDir[dir](i);
      }

    }
    t2 = clock();
    tm1 += (double) (t2-t1)/CLOCKS_PER_SEC;

    // Check to see if we need to retake time step
    if (cfla > cflm)
      return Lucee::UpdaterStatus(false, dt*cfl/cfla);

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
      double vCoord[3];
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

          q.setPtr(qPtr, idx);
          q.setPtr(qPtrl, idxl);
          EM.setPtr(emPtr, idx);
          EM.setPtr(emPtrl, idxl);

          // Copy data to Eigen vectors
          for (int i = 0; i < nlocal; i++)
          {
            rightData(i) = qPtr[i];
            leftData(i) = qPtrl[i];
          }

          rightDataAtQuad.noalias() = surfLowerQuad[dir].interpMat*rightData;
          leftDataAtQuad.noalias() = surfUpperQuad[dir].interpMat*leftData;
          calcFlux(dir, phaseNodeCoords, emPtr, flux);
          alphaRight.noalias() = surfLowerQuad[dir].interpMat*flux;
          calcFlux(dir, phaseNodeCoords, emPtrl, flux);
          alphaLeft.noalias() = surfUpperQuad[dir].interpMat*flux;

          // Compute numerical flux
          computeNumericalFlux(alphaLeft, alphaRight, leftDataAtQuad, rightDataAtQuad, numericalFluxAtQuad);

          qNew.setPtr(qNewPtr, idx);
          qNew.setPtr(qNewPtrl, idxl);

          resultVectorDirLeft[dir].noalias() = bigStoredUpperSurfMatrices[dir]*(numericalFluxAtQuad.cwiseProduct(surfUpperQuad[dir].weights));
          resultVectorDir[dir].noalias() = bigStoredLowerSurfMatrices[dir]*(numericalFluxAtQuad.cwiseProduct(surfLowerQuad[dir].weights));

          for (int i = 0; i < nlocal; i++)
          {
            qNewPtrl[i] -= resultVectorDirLeft[dir](i);
            qNewPtr[i] += resultVectorDir[dir](i);
          }
        }
      }
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
  EigenNodalVlasovUpdater<CDIM,VDIM>::computeNumericalFlux(const Eigen::VectorXd& alphaLeft,
    const Eigen::VectorXd& alphaRight, const Eigen::VectorXd& leftValsAtQuad, 
    const Eigen::VectorXd& rightValsAtQuad, Eigen::VectorXd& numericalFluxAtQuad)
  {
    double maxs;
    // Loop through all quadrature points
    for (int quadIndex = 0; quadIndex < numericalFluxAtQuad.size(); quadIndex++)
    {
      maxs = std::max(std::fabs(alphaLeft(quadIndex)),std::fabs(alphaRight(quadIndex)));
      numericalFluxAtQuad(quadIndex) = 0.5*(alphaLeft(quadIndex)*leftValsAtQuad(quadIndex) + alphaRight(quadIndex)*rightValsAtQuad(quadIndex)) 
        - 0.5*maxs*(rightValsAtQuad(quadIndex) - leftValsAtQuad(quadIndex));

    }
  }

  template <unsigned CDIM, unsigned VDIM>
  void
  EigenNodalVlasovUpdater<CDIM,VDIM>::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
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
