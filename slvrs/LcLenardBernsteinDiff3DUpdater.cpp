/**
 * @file	LcLenardBernsteinDiff3DUpdater.cpp
 *
 * @brief	Updater to evaluate the diffusion term in the L-B collision operator for 1D2V problems.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLenardBernsteinDiff3DUpdater.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  using namespace Eigen;
  const char *LenardBernsteinDiff3DUpdater::id = "LenardBernsteinDiff3DUpdater";

  LenardBernsteinDiff3DUpdater::LenardBernsteinDiff3DUpdater()
  {
  }

  LenardBernsteinDiff3DUpdater::~LenardBernsteinDiff3DUpdater()
  {
  }

  void
  LenardBernsteinDiff3DUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis");
    else
      throw Lucee::Except("LenardBernsteinDiff3DUpdater::readInput: Must specify element to use using 'basis'");
 
    // CFL number to control time-step
    cfl = tbl.getNumber("cfl"); // CFL number
    // should only increments be computed?
    onlyIncrement = false;
    if (tbl.hasBool("onlyIncrement"))
      onlyIncrement = tbl.getBool("onlyIncrement");

    if (tbl.hasNumber("speciesMass"))
      speciesMass = tbl.getNumber("speciesMass");
    else
      throw Lucee::Except("LenardBernsteinDiff3DUpdater::readInput: Must specify speciesMass");

    if (tbl.hasNumber("B0"))
      B0 = tbl.getNumber("B0");
    else
      throw Lucee::Except("LenardBernsteinDiff3DUpdater::readInput: Must specify magnetic field B0");

    if (tbl.hasFunction("alpha"))
      fnRef = tbl.getFunctionRef("alpha");
    else
      throw Lucee::Except("LenardBernsteinDiff3DUpdater::readInput: Must supply a collision frequency function as alpha.");
  
  }

  void
  LenardBernsteinDiff3DUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();
    
// get hold of grid
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();
    // local region to update
    Lucee::Region<3, int> localRgn = grid.getLocalRegion();

    unsigned nlocal = nodalBasis->getNumNodes();
    // allocate space for matrices
    std::vector<Lucee::Matrix<double> > iMatLucee;
    std::vector<Lucee::Matrix<double> > lowerMatLucee;
    std::vector<Lucee::Matrix<double> > upperMatLucee;
    iMatLucee.resize(3);
    lowerMatLucee.resize(3);
    upperMatLucee.resize(3);

    iMat.resize(3);
    lowerMat.resize(3);
    upperMat.resize(3);

    for (int dir = 0; dir < 3; dir++)
    {
      iMatLucee[dir] = Lucee::Matrix<double>(nlocal, nlocal);
      lowerMatLucee[dir] = Lucee::Matrix<double>(nlocal, nlocal);
      upperMatLucee[dir] = Lucee::Matrix<double>(nlocal, nlocal);

      iMat[dir] = Eigen::MatrixXd(nlocal, nlocal);
      lowerMat[dir] = Eigen::MatrixXd(nlocal, nlocal);
      upperMat[dir] = Eigen::MatrixXd(nlocal, nlocal);
    }

    Lucee::RowMajorSequencer<3> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[3];
    seq.fillWithIndex(idx);
    nodalBasis->setIndex(idx);

    // get various matrices needed
    nodalBasis->getDiffusionMatrices(iMatLucee, lowerMatLucee, upperMatLucee);
    for (int dir = 0; dir < 3; dir ++)
    {
      copyLuceeToEigen(iMatLucee[dir], iMat[dir]);
      copyLuceeToEigen(lowerMatLucee[dir], lowerMat[dir]);
      copyLuceeToEigen(upperMatLucee[dir], upperMat[dir]);
    }

    // pre-multiply each of the matrices by inverse matrix
    Lucee::Matrix<double> massMatrixLucee(nlocal, nlocal);
    // NOTE: mass matrix is fetched repeatedly as the solve() method
    // destroys it during the inversion process
    nodalBasis->getMassMatrix(massMatrixLucee);
    Eigen::MatrixXd massMatrix(nlocal, nlocal);
    copyLuceeToEigen(massMatrixLucee, massMatrix);
    Eigen::MatrixXd massMatrixInv = massMatrix.inverse();

    for (int dir = 0; dir < 3; dir++)
    {
      iMat[dir] = massMatrixInv*iMat[dir];
      lowerMat[dir] = massMatrixInv*lowerMat[dir];
      upperMat[dir] = massMatrixInv*upperMat[dir];
    }

    int numSurfQuadNodes = nodalBasis->getNumSurfGaussNodes();
    Lucee::Matrix<double> interpSurfMatrixLowerLucee(numSurfQuadNodes, nlocal);
    Lucee::Matrix<double> gaussSurfOrdinatesLucee(numSurfQuadNodes, 3);
    gaussSurfWeights = std::vector<double>(numSurfQuadNodes);
    gaussSurfOrdinates = Eigen::MatrixXd(numSurfQuadNodes, 3);
    std::vector<int> lowerSurfNodeNums(nodalBasis->getNumSurfLowerNodes(1));
    Eigen::MatrixXd interpSurfMatrixLower(numSurfQuadNodes, nlocal);
    
    // Get the interpolation matrix for the lower surface quadrature points.
    nodalBasis->getSurfLowerGaussQuadData(1, interpSurfMatrixLowerLucee, gaussSurfOrdinatesLucee,
      gaussSurfWeights);
    // Get the nodes on the lower surface. Use this with interpSurfMatrixLucee.
    nodalBasis->getSurfLowerNodeNums(1, lowerSurfNodeNums);
    // Matrix to compute v_t(x)^2 at quadrature locations. Only works in 2-D.
    surfNodeInterpMatrix = Eigen::MatrixXd(numSurfQuadNodes, lowerSurfNodeNums.size());

    copyLuceeToEigen(interpSurfMatrixLowerLucee, interpSurfMatrixLower);

    // Take interpSurfMatrixLower and create a lower dimension interpolation matrix
    for (int nodeIndex = 0; nodeIndex < numSurfQuadNodes; nodeIndex++)
    {
      // At each quadrature node, copy basis function evaluations for
      // those basis functions associated with the nodes on the lower surface
      for (int basisIndex = 0; basisIndex < lowerSurfNodeNums.size(); basisIndex++)
      {
        // Need order of elements in lowerSurfNodeNums to match up with the
        // 1-D data in u(x) for this to work. Not sure how to enforce, but maybe
        // it will just happen to be that way.
        surfNodeInterpMatrix(nodeIndex, basisIndex) = interpSurfMatrixLower(nodeIndex, 
          lowerSurfNodeNums[basisIndex]);
      }
    }

    // Additional matrices
    iMatDiffusionTimesMu = Eigen::MatrixXd(nlocal, nlocal);
    lowerMatDiffusionTimesMu = Eigen::MatrixXd(nlocal, nlocal);
    upperMatDiffusionTimesMu = Eigen::MatrixXd(nlocal, nlocal);

    double dq[3];
    for (int i = 0; i < 3; i++)
      dq[i] = grid.getDx(i);

    // First order serendipity element
    if (nlocal == 8)
    {
      iMatDiffusionTimesMu(0,0)=(5*dq[0]*dq[1])/72;
      iMatDiffusionTimesMu(0,1)=(5*dq[0]*dq[1])/144;
      iMatDiffusionTimesMu(0,2)=(5*dq[0]*dq[1])/144;
      iMatDiffusionTimesMu(0,3)=(5*dq[0]*dq[1])/288;
      iMatDiffusionTimesMu(0,4)=0;
      iMatDiffusionTimesMu(0,5)=0;
      iMatDiffusionTimesMu(0,6)=0;
      iMatDiffusionTimesMu(0,7)=0;
      iMatDiffusionTimesMu(1,0)=(5*dq[0]*dq[1])/144;
      iMatDiffusionTimesMu(1,1)=(5*dq[0]*dq[1])/72;
      iMatDiffusionTimesMu(1,2)=(5*dq[0]*dq[1])/288;
      iMatDiffusionTimesMu(1,3)=(5*dq[0]*dq[1])/144;
      iMatDiffusionTimesMu(1,4)=0;
      iMatDiffusionTimesMu(1,5)=0;
      iMatDiffusionTimesMu(1,6)=0;
      iMatDiffusionTimesMu(1,7)=0;
      iMatDiffusionTimesMu(2,0)=(5*dq[0]*dq[1])/144;
      iMatDiffusionTimesMu(2,1)=(5*dq[0]*dq[1])/288;
      iMatDiffusionTimesMu(2,2)=(5*dq[0]*dq[1])/72;
      iMatDiffusionTimesMu(2,3)=(5*dq[0]*dq[1])/144;
      iMatDiffusionTimesMu(2,4)=0;
      iMatDiffusionTimesMu(2,5)=0;
      iMatDiffusionTimesMu(2,6)=0;
      iMatDiffusionTimesMu(2,7)=0;
      iMatDiffusionTimesMu(3,0)=(5*dq[0]*dq[1])/288;
      iMatDiffusionTimesMu(3,1)=(5*dq[0]*dq[1])/144;
      iMatDiffusionTimesMu(3,2)=(5*dq[0]*dq[1])/144;
      iMatDiffusionTimesMu(3,3)=(5*dq[0]*dq[1])/72;
      iMatDiffusionTimesMu(3,4)=0;
      iMatDiffusionTimesMu(3,5)=0;
      iMatDiffusionTimesMu(3,6)=0;
      iMatDiffusionTimesMu(3,7)=0;
      iMatDiffusionTimesMu(4,0)=0;
      iMatDiffusionTimesMu(4,1)=0;
      iMatDiffusionTimesMu(4,2)=0;
      iMatDiffusionTimesMu(4,3)=0;
      iMatDiffusionTimesMu(4,4)=(-5*dq[0]*dq[1])/72;
      iMatDiffusionTimesMu(4,5)=(-5*dq[0]*dq[1])/144;
      iMatDiffusionTimesMu(4,6)=(-5*dq[0]*dq[1])/144;
      iMatDiffusionTimesMu(4,7)=(-5*dq[0]*dq[1])/288;
      iMatDiffusionTimesMu(5,0)=0;
      iMatDiffusionTimesMu(5,1)=0;
      iMatDiffusionTimesMu(5,2)=0;
      iMatDiffusionTimesMu(5,3)=0;
      iMatDiffusionTimesMu(5,4)=(-5*dq[0]*dq[1])/144;
      iMatDiffusionTimesMu(5,5)=(-5*dq[0]*dq[1])/72;
      iMatDiffusionTimesMu(5,6)=(-5*dq[0]*dq[1])/288;
      iMatDiffusionTimesMu(5,7)=(-5*dq[0]*dq[1])/144;
      iMatDiffusionTimesMu(6,0)=0;
      iMatDiffusionTimesMu(6,1)=0;
      iMatDiffusionTimesMu(6,2)=0;
      iMatDiffusionTimesMu(6,3)=0;
      iMatDiffusionTimesMu(6,4)=(-5*dq[0]*dq[1])/144;
      iMatDiffusionTimesMu(6,5)=(-5*dq[0]*dq[1])/288;
      iMatDiffusionTimesMu(6,6)=(-5*dq[0]*dq[1])/72;
      iMatDiffusionTimesMu(6,7)=(-5*dq[0]*dq[1])/144;
      iMatDiffusionTimesMu(7,0)=0;
      iMatDiffusionTimesMu(7,1)=0;
      iMatDiffusionTimesMu(7,2)=0;
      iMatDiffusionTimesMu(7,3)=0;
      iMatDiffusionTimesMu(7,4)=(-5*dq[0]*dq[1])/288;
      iMatDiffusionTimesMu(7,5)=(-5*dq[0]*dq[1])/144;
      iMatDiffusionTimesMu(7,6)=(-5*dq[0]*dq[1])/144;
      iMatDiffusionTimesMu(7,7)=(-5*dq[0]*dq[1])/72;
      lowerMatDiffusionTimesMu(0,0)=(-5*dq[0]*dq[1])/216;
      lowerMatDiffusionTimesMu(0,1)=(-5*dq[0]*dq[1])/432;
      lowerMatDiffusionTimesMu(0,2)=(-5*dq[0]*dq[1])/432;
      lowerMatDiffusionTimesMu(0,3)=(-5*dq[0]*dq[1])/864;
      lowerMatDiffusionTimesMu(0,4)=(-2*dq[0]*dq[1])/27;
      lowerMatDiffusionTimesMu(0,5)=-(dq[0]*dq[1])/27;
      lowerMatDiffusionTimesMu(0,6)=-(dq[0]*dq[1])/27;
      lowerMatDiffusionTimesMu(0,7)=-(dq[0]*dq[1])/54;
      lowerMatDiffusionTimesMu(1,0)=(-5*dq[0]*dq[1])/432;
      lowerMatDiffusionTimesMu(1,1)=(-5*dq[0]*dq[1])/216;
      lowerMatDiffusionTimesMu(1,2)=(-5*dq[0]*dq[1])/864;
      lowerMatDiffusionTimesMu(1,3)=(-5*dq[0]*dq[1])/432;
      lowerMatDiffusionTimesMu(1,4)=-(dq[0]*dq[1])/27;
      lowerMatDiffusionTimesMu(1,5)=(-2*dq[0]*dq[1])/27;
      lowerMatDiffusionTimesMu(1,6)=-(dq[0]*dq[1])/54;
      lowerMatDiffusionTimesMu(1,7)=-(dq[0]*dq[1])/27;
      lowerMatDiffusionTimesMu(2,0)=(-5*dq[0]*dq[1])/432;
      lowerMatDiffusionTimesMu(2,1)=(-5*dq[0]*dq[1])/864;
      lowerMatDiffusionTimesMu(2,2)=(-5*dq[0]*dq[1])/216;
      lowerMatDiffusionTimesMu(2,3)=(-5*dq[0]*dq[1])/432;
      lowerMatDiffusionTimesMu(2,4)=-(dq[0]*dq[1])/27;
      lowerMatDiffusionTimesMu(2,5)=-(dq[0]*dq[1])/54;
      lowerMatDiffusionTimesMu(2,6)=(-2*dq[0]*dq[1])/27;
      lowerMatDiffusionTimesMu(2,7)=-(dq[0]*dq[1])/27;
      lowerMatDiffusionTimesMu(3,0)=(-5*dq[0]*dq[1])/864;
      lowerMatDiffusionTimesMu(3,1)=(-5*dq[0]*dq[1])/432;
      lowerMatDiffusionTimesMu(3,2)=(-5*dq[0]*dq[1])/432;
      lowerMatDiffusionTimesMu(3,3)=(-5*dq[0]*dq[1])/216;
      lowerMatDiffusionTimesMu(3,4)=-(dq[0]*dq[1])/54;
      lowerMatDiffusionTimesMu(3,5)=-(dq[0]*dq[1])/27;
      lowerMatDiffusionTimesMu(3,6)=-(dq[0]*dq[1])/27;
      lowerMatDiffusionTimesMu(3,7)=(-2*dq[0]*dq[1])/27;
      lowerMatDiffusionTimesMu(4,0)=-(dq[0]*dq[1])/216;
      lowerMatDiffusionTimesMu(4,1)=-(dq[0]*dq[1])/432;
      lowerMatDiffusionTimesMu(4,2)=-(dq[0]*dq[1])/432;
      lowerMatDiffusionTimesMu(4,3)=-(dq[0]*dq[1])/864;
      lowerMatDiffusionTimesMu(4,4)=(-5*dq[0]*dq[1])/216;
      lowerMatDiffusionTimesMu(4,5)=(-5*dq[0]*dq[1])/432;
      lowerMatDiffusionTimesMu(4,6)=(-5*dq[0]*dq[1])/432;
      lowerMatDiffusionTimesMu(4,7)=(-5*dq[0]*dq[1])/864;
      lowerMatDiffusionTimesMu(5,0)=-(dq[0]*dq[1])/432;
      lowerMatDiffusionTimesMu(5,1)=-(dq[0]*dq[1])/216;
      lowerMatDiffusionTimesMu(5,2)=-(dq[0]*dq[1])/864;
      lowerMatDiffusionTimesMu(5,3)=-(dq[0]*dq[1])/432;
      lowerMatDiffusionTimesMu(5,4)=(-5*dq[0]*dq[1])/432;
      lowerMatDiffusionTimesMu(5,5)=(-5*dq[0]*dq[1])/216;
      lowerMatDiffusionTimesMu(5,6)=(-5*dq[0]*dq[1])/864;
      lowerMatDiffusionTimesMu(5,7)=(-5*dq[0]*dq[1])/432;
      lowerMatDiffusionTimesMu(6,0)=-(dq[0]*dq[1])/432;
      lowerMatDiffusionTimesMu(6,1)=-(dq[0]*dq[1])/864;
      lowerMatDiffusionTimesMu(6,2)=-(dq[0]*dq[1])/216;
      lowerMatDiffusionTimesMu(6,3)=-(dq[0]*dq[1])/432;
      lowerMatDiffusionTimesMu(6,4)=(-5*dq[0]*dq[1])/432;
      lowerMatDiffusionTimesMu(6,5)=(-5*dq[0]*dq[1])/864;
      lowerMatDiffusionTimesMu(6,6)=(-5*dq[0]*dq[1])/216;
      lowerMatDiffusionTimesMu(6,7)=(-5*dq[0]*dq[1])/432;
      lowerMatDiffusionTimesMu(7,0)=-(dq[0]*dq[1])/864;
      lowerMatDiffusionTimesMu(7,1)=-(dq[0]*dq[1])/432;
      lowerMatDiffusionTimesMu(7,2)=-(dq[0]*dq[1])/432;
      lowerMatDiffusionTimesMu(7,3)=-(dq[0]*dq[1])/216;
      lowerMatDiffusionTimesMu(7,4)=(-5*dq[0]*dq[1])/864;
      lowerMatDiffusionTimesMu(7,5)=(-5*dq[0]*dq[1])/432;
      lowerMatDiffusionTimesMu(7,6)=(-5*dq[0]*dq[1])/432;
      lowerMatDiffusionTimesMu(7,7)=(-5*dq[0]*dq[1])/216;
      upperMatDiffusionTimesMu(0,0)=(5*dq[0]*dq[1])/216;
      upperMatDiffusionTimesMu(0,1)=(5*dq[0]*dq[1])/432;
      upperMatDiffusionTimesMu(0,2)=(5*dq[0]*dq[1])/432;
      upperMatDiffusionTimesMu(0,3)=(5*dq[0]*dq[1])/864;
      upperMatDiffusionTimesMu(0,4)=(dq[0]*dq[1])/216;
      upperMatDiffusionTimesMu(0,5)=(dq[0]*dq[1])/432;
      upperMatDiffusionTimesMu(0,6)=(dq[0]*dq[1])/432;
      upperMatDiffusionTimesMu(0,7)=(dq[0]*dq[1])/864;
      upperMatDiffusionTimesMu(1,0)=(5*dq[0]*dq[1])/432;
      upperMatDiffusionTimesMu(1,1)=(5*dq[0]*dq[1])/216;
      upperMatDiffusionTimesMu(1,2)=(5*dq[0]*dq[1])/864;
      upperMatDiffusionTimesMu(1,3)=(5*dq[0]*dq[1])/432;
      upperMatDiffusionTimesMu(1,4)=(dq[0]*dq[1])/432;
      upperMatDiffusionTimesMu(1,5)=(dq[0]*dq[1])/216;
      upperMatDiffusionTimesMu(1,6)=(dq[0]*dq[1])/864;
      upperMatDiffusionTimesMu(1,7)=(dq[0]*dq[1])/432;
      upperMatDiffusionTimesMu(2,0)=(5*dq[0]*dq[1])/432;
      upperMatDiffusionTimesMu(2,1)=(5*dq[0]*dq[1])/864;
      upperMatDiffusionTimesMu(2,2)=(5*dq[0]*dq[1])/216;
      upperMatDiffusionTimesMu(2,3)=(5*dq[0]*dq[1])/432;
      upperMatDiffusionTimesMu(2,4)=(dq[0]*dq[1])/432;
      upperMatDiffusionTimesMu(2,5)=(dq[0]*dq[1])/864;
      upperMatDiffusionTimesMu(2,6)=(dq[0]*dq[1])/216;
      upperMatDiffusionTimesMu(2,7)=(dq[0]*dq[1])/432;
      upperMatDiffusionTimesMu(3,0)=(5*dq[0]*dq[1])/864;
      upperMatDiffusionTimesMu(3,1)=(5*dq[0]*dq[1])/432;
      upperMatDiffusionTimesMu(3,2)=(5*dq[0]*dq[1])/432;
      upperMatDiffusionTimesMu(3,3)=(5*dq[0]*dq[1])/216;
      upperMatDiffusionTimesMu(3,4)=(dq[0]*dq[1])/864;
      upperMatDiffusionTimesMu(3,5)=(dq[0]*dq[1])/432;
      upperMatDiffusionTimesMu(3,6)=(dq[0]*dq[1])/432;
      upperMatDiffusionTimesMu(3,7)=(dq[0]*dq[1])/216;
      upperMatDiffusionTimesMu(4,0)=(2*dq[0]*dq[1])/27;
      upperMatDiffusionTimesMu(4,1)=(dq[0]*dq[1])/27;
      upperMatDiffusionTimesMu(4,2)=(dq[0]*dq[1])/27;
      upperMatDiffusionTimesMu(4,3)=(dq[0]*dq[1])/54;
      upperMatDiffusionTimesMu(4,4)=(5*dq[0]*dq[1])/216;
      upperMatDiffusionTimesMu(4,5)=(5*dq[0]*dq[1])/432;
      upperMatDiffusionTimesMu(4,6)=(5*dq[0]*dq[1])/432;
      upperMatDiffusionTimesMu(4,7)=(5*dq[0]*dq[1])/864;
      upperMatDiffusionTimesMu(5,0)=(dq[0]*dq[1])/27;
      upperMatDiffusionTimesMu(5,1)=(2*dq[0]*dq[1])/27;
      upperMatDiffusionTimesMu(5,2)=(dq[0]*dq[1])/54;
      upperMatDiffusionTimesMu(5,3)=(dq[0]*dq[1])/27;
      upperMatDiffusionTimesMu(5,4)=(5*dq[0]*dq[1])/432;
      upperMatDiffusionTimesMu(5,5)=(5*dq[0]*dq[1])/216;
      upperMatDiffusionTimesMu(5,6)=(5*dq[0]*dq[1])/864;
      upperMatDiffusionTimesMu(5,7)=(5*dq[0]*dq[1])/432;
      upperMatDiffusionTimesMu(6,0)=(dq[0]*dq[1])/27;
      upperMatDiffusionTimesMu(6,1)=(dq[0]*dq[1])/54;
      upperMatDiffusionTimesMu(6,2)=(2*dq[0]*dq[1])/27;
      upperMatDiffusionTimesMu(6,3)=(dq[0]*dq[1])/27;
      upperMatDiffusionTimesMu(6,4)=(5*dq[0]*dq[1])/432;
      upperMatDiffusionTimesMu(6,5)=(5*dq[0]*dq[1])/864;
      upperMatDiffusionTimesMu(6,6)=(5*dq[0]*dq[1])/216;
      upperMatDiffusionTimesMu(6,7)=(5*dq[0]*dq[1])/432;
      upperMatDiffusionTimesMu(7,0)=(dq[0]*dq[1])/54;
      upperMatDiffusionTimesMu(7,1)=(dq[0]*dq[1])/27;
      upperMatDiffusionTimesMu(7,2)=(dq[0]*dq[1])/27;
      upperMatDiffusionTimesMu(7,3)=(2*dq[0]*dq[1])/27;
      upperMatDiffusionTimesMu(7,4)=(5*dq[0]*dq[1])/864;
      upperMatDiffusionTimesMu(7,5)=(5*dq[0]*dq[1])/432;
      upperMatDiffusionTimesMu(7,6)=(5*dq[0]*dq[1])/432;
      upperMatDiffusionTimesMu(7,7)=(5*dq[0]*dq[1])/216;
    }
    else if (nlocal == 20)
    {
      iMatDiffusionTimesMu(0,0)=(11*dq[0]*dq[1])/240;
      iMatDiffusionTimesMu(0,1)=(-323*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(0,2)=(239*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(0,3)=(-323*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(0,4)=(-653*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(0,5)=(239*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(0,6)=(-653*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(0,7)=(61*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(0,8)=0;
      iMatDiffusionTimesMu(0,9)=(-89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(0,10)=(-89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(0,11)=(-89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(0,12)=0;
      iMatDiffusionTimesMu(0,13)=(25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(0,14)=-(dq[0]*dq[1])/480;
      iMatDiffusionTimesMu(0,15)=(25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(0,16)=(25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(0,17)=-(dq[0]*dq[1])/480;
      iMatDiffusionTimesMu(0,18)=(25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(0,19)=-(dq[0]*dq[1])/480;
      iMatDiffusionTimesMu(1,0)=(-287*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(1,1)=(11*dq[0]*dq[1])/45;
      iMatDiffusionTimesMu(1,2)=(-287*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(1,3)=(11*dq[0]*dq[1])/72;
      iMatDiffusionTimesMu(1,4)=(11*dq[0]*dq[1])/72;
      iMatDiffusionTimesMu(1,5)=(-617*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(1,6)=(11*dq[0]*dq[1])/90;
      iMatDiffusionTimesMu(1,7)=(-617*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(1,8)=(89*dq[0]*dq[1])/2160;
      iMatDiffusionTimesMu(1,9)=(89*dq[0]*dq[1])/2160;
      iMatDiffusionTimesMu(1,10)=(89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(1,11)=(89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(1,12)=(-89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(1,13)=0;
      iMatDiffusionTimesMu(1,14)=(-89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(1,15)=0;
      iMatDiffusionTimesMu(1,16)=0;
      iMatDiffusionTimesMu(1,17)=(-89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(1,18)=0;
      iMatDiffusionTimesMu(1,19)=(-89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(2,0)=(239*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(2,1)=(-323*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(2,2)=(11*dq[0]*dq[1])/240;
      iMatDiffusionTimesMu(2,3)=(-653*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(2,4)=(-323*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(2,5)=(61*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(2,6)=(-653*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(2,7)=(239*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(2,8)=(-89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(2,9)=0;
      iMatDiffusionTimesMu(2,10)=(-89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(2,11)=(-89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(2,12)=-(dq[0]*dq[1])/480;
      iMatDiffusionTimesMu(2,13)=(25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(2,14)=0;
      iMatDiffusionTimesMu(2,15)=(25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(2,16)=(25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(2,17)=-(dq[0]*dq[1])/480;
      iMatDiffusionTimesMu(2,18)=(25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(2,19)=-(dq[0]*dq[1])/480;
      iMatDiffusionTimesMu(3,0)=(-287*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(3,1)=(11*dq[0]*dq[1])/72;
      iMatDiffusionTimesMu(3,2)=(-617*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(3,3)=(11*dq[0]*dq[1])/45;
      iMatDiffusionTimesMu(3,4)=(11*dq[0]*dq[1])/90;
      iMatDiffusionTimesMu(3,5)=(-287*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(3,6)=(11*dq[0]*dq[1])/72;
      iMatDiffusionTimesMu(3,7)=(-617*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(3,8)=(89*dq[0]*dq[1])/2160;
      iMatDiffusionTimesMu(3,9)=(89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(3,10)=(89*dq[0]*dq[1])/2160;
      iMatDiffusionTimesMu(3,11)=(89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(3,12)=(-89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(3,13)=0;
      iMatDiffusionTimesMu(3,14)=(-89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(3,15)=0;
      iMatDiffusionTimesMu(3,16)=0;
      iMatDiffusionTimesMu(3,17)=(-89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(3,18)=0;
      iMatDiffusionTimesMu(3,19)=(-89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(4,0)=(-617*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(4,1)=(11*dq[0]*dq[1])/72;
      iMatDiffusionTimesMu(4,2)=(-287*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(4,3)=(11*dq[0]*dq[1])/90;
      iMatDiffusionTimesMu(4,4)=(11*dq[0]*dq[1])/45;
      iMatDiffusionTimesMu(4,5)=(-617*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(4,6)=(11*dq[0]*dq[1])/72;
      iMatDiffusionTimesMu(4,7)=(-287*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(4,8)=(89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(4,9)=(89*dq[0]*dq[1])/2160;
      iMatDiffusionTimesMu(4,10)=(89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(4,11)=(89*dq[0]*dq[1])/2160;
      iMatDiffusionTimesMu(4,12)=(-89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(4,13)=0;
      iMatDiffusionTimesMu(4,14)=(-89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(4,15)=0;
      iMatDiffusionTimesMu(4,16)=0;
      iMatDiffusionTimesMu(4,17)=(-89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(4,18)=0;
      iMatDiffusionTimesMu(4,19)=(-89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(5,0)=(239*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(5,1)=(-653*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(5,2)=(61*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(5,3)=(-323*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(5,4)=(-653*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(5,5)=(11*dq[0]*dq[1])/240;
      iMatDiffusionTimesMu(5,6)=(-323*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(5,7)=(239*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(5,8)=(-89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(5,9)=(-89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(5,10)=0;
      iMatDiffusionTimesMu(5,11)=(-89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(5,12)=-(dq[0]*dq[1])/480;
      iMatDiffusionTimesMu(5,13)=(25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(5,14)=-(dq[0]*dq[1])/480;
      iMatDiffusionTimesMu(5,15)=(25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(5,16)=(25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(5,17)=0;
      iMatDiffusionTimesMu(5,18)=(25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(5,19)=-(dq[0]*dq[1])/480;
      iMatDiffusionTimesMu(6,0)=(-617*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(6,1)=(11*dq[0]*dq[1])/90;
      iMatDiffusionTimesMu(6,2)=(-617*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(6,3)=(11*dq[0]*dq[1])/72;
      iMatDiffusionTimesMu(6,4)=(11*dq[0]*dq[1])/72;
      iMatDiffusionTimesMu(6,5)=(-287*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(6,6)=(11*dq[0]*dq[1])/45;
      iMatDiffusionTimesMu(6,7)=(-287*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(6,8)=(89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(6,9)=(89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(6,10)=(89*dq[0]*dq[1])/2160;
      iMatDiffusionTimesMu(6,11)=(89*dq[0]*dq[1])/2160;
      iMatDiffusionTimesMu(6,12)=(-89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(6,13)=0;
      iMatDiffusionTimesMu(6,14)=(-89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(6,15)=0;
      iMatDiffusionTimesMu(6,16)=0;
      iMatDiffusionTimesMu(6,17)=(-89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(6,18)=0;
      iMatDiffusionTimesMu(6,19)=(-89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(7,0)=(61*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(7,1)=(-653*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(7,2)=(239*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(7,3)=(-653*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(7,4)=(-323*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(7,5)=(239*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(7,6)=(-323*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(7,7)=(11*dq[0]*dq[1])/240;
      iMatDiffusionTimesMu(7,8)=(-89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(7,9)=(-89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(7,10)=(-89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(7,11)=0;
      iMatDiffusionTimesMu(7,12)=-(dq[0]*dq[1])/480;
      iMatDiffusionTimesMu(7,13)=(25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(7,14)=-(dq[0]*dq[1])/480;
      iMatDiffusionTimesMu(7,15)=(25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(7,16)=(25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(7,17)=-(dq[0]*dq[1])/480;
      iMatDiffusionTimesMu(7,18)=(25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(7,19)=0;
      iMatDiffusionTimesMu(8,0)=0;
      iMatDiffusionTimesMu(8,1)=(25*dq[0]*dq[1])/432;
      iMatDiffusionTimesMu(8,2)=(-25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(8,3)=(25*dq[0]*dq[1])/432;
      iMatDiffusionTimesMu(8,4)=(25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(8,5)=(-25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(8,6)=(25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(8,7)=(-25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(8,8)=0;
      iMatDiffusionTimesMu(8,9)=0;
      iMatDiffusionTimesMu(8,10)=0;
      iMatDiffusionTimesMu(8,11)=0;
      iMatDiffusionTimesMu(8,12)=0;
      iMatDiffusionTimesMu(8,13)=(-25*dq[0]*dq[1])/432;
      iMatDiffusionTimesMu(8,14)=(25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(8,15)=(-25*dq[0]*dq[1])/432;
      iMatDiffusionTimesMu(8,16)=(-25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(8,17)=(25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(8,18)=(-25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(8,19)=(25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(9,0)=(-25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(9,1)=(25*dq[0]*dq[1])/432;
      iMatDiffusionTimesMu(9,2)=0;
      iMatDiffusionTimesMu(9,3)=(25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(9,4)=(25*dq[0]*dq[1])/432;
      iMatDiffusionTimesMu(9,5)=(-25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(9,6)=(25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(9,7)=(-25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(9,8)=0;
      iMatDiffusionTimesMu(9,9)=0;
      iMatDiffusionTimesMu(9,10)=0;
      iMatDiffusionTimesMu(9,11)=0;
      iMatDiffusionTimesMu(9,12)=(25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(9,13)=(-25*dq[0]*dq[1])/432;
      iMatDiffusionTimesMu(9,14)=0;
      iMatDiffusionTimesMu(9,15)=(-25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(9,16)=(-25*dq[0]*dq[1])/432;
      iMatDiffusionTimesMu(9,17)=(25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(9,18)=(-25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(9,19)=(25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(10,0)=(-25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(10,1)=(25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(10,2)=(-25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(10,3)=(25*dq[0]*dq[1])/432;
      iMatDiffusionTimesMu(10,4)=(25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(10,5)=0;
      iMatDiffusionTimesMu(10,6)=(25*dq[0]*dq[1])/432;
      iMatDiffusionTimesMu(10,7)=(-25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(10,8)=0;
      iMatDiffusionTimesMu(10,9)=0;
      iMatDiffusionTimesMu(10,10)=0;
      iMatDiffusionTimesMu(10,11)=0;
      iMatDiffusionTimesMu(10,12)=(25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(10,13)=(-25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(10,14)=(25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(10,15)=(-25*dq[0]*dq[1])/432;
      iMatDiffusionTimesMu(10,16)=(-25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(10,17)=0;
      iMatDiffusionTimesMu(10,18)=(-25*dq[0]*dq[1])/432;
      iMatDiffusionTimesMu(10,19)=(25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(11,0)=(-25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(11,1)=(25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(11,2)=(-25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(11,3)=(25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(11,4)=(25*dq[0]*dq[1])/432;
      iMatDiffusionTimesMu(11,5)=(-25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(11,6)=(25*dq[0]*dq[1])/432;
      iMatDiffusionTimesMu(11,7)=0;
      iMatDiffusionTimesMu(11,8)=0;
      iMatDiffusionTimesMu(11,9)=0;
      iMatDiffusionTimesMu(11,10)=0;
      iMatDiffusionTimesMu(11,11)=0;
      iMatDiffusionTimesMu(11,12)=(25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(11,13)=(-25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(11,14)=(25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(11,15)=(-25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(11,16)=(-25*dq[0]*dq[1])/432;
      iMatDiffusionTimesMu(11,17)=(25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(11,18)=(-25*dq[0]*dq[1])/432;
      iMatDiffusionTimesMu(11,19)=0;
      iMatDiffusionTimesMu(12,0)=0;
      iMatDiffusionTimesMu(12,1)=(-25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(12,2)=(dq[0]*dq[1])/480;
      iMatDiffusionTimesMu(12,3)=(-25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(12,4)=(-25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(12,5)=(dq[0]*dq[1])/480;
      iMatDiffusionTimesMu(12,6)=(-25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(12,7)=(dq[0]*dq[1])/480;
      iMatDiffusionTimesMu(12,8)=0;
      iMatDiffusionTimesMu(12,9)=(89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(12,10)=(89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(12,11)=(89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(12,12)=(-11*dq[0]*dq[1])/240;
      iMatDiffusionTimesMu(12,13)=(323*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(12,14)=(-239*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(12,15)=(323*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(12,16)=(653*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(12,17)=(-239*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(12,18)=(653*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(12,19)=(-61*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(13,0)=(89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(13,1)=0;
      iMatDiffusionTimesMu(13,2)=(89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(13,3)=0;
      iMatDiffusionTimesMu(13,4)=0;
      iMatDiffusionTimesMu(13,5)=(89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(13,6)=0;
      iMatDiffusionTimesMu(13,7)=(89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(13,8)=(-89*dq[0]*dq[1])/2160;
      iMatDiffusionTimesMu(13,9)=(-89*dq[0]*dq[1])/2160;
      iMatDiffusionTimesMu(13,10)=(-89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(13,11)=(-89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(13,12)=(287*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(13,13)=(-11*dq[0]*dq[1])/45;
      iMatDiffusionTimesMu(13,14)=(287*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(13,15)=(-11*dq[0]*dq[1])/72;
      iMatDiffusionTimesMu(13,16)=(-11*dq[0]*dq[1])/72;
      iMatDiffusionTimesMu(13,17)=(617*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(13,18)=(-11*dq[0]*dq[1])/90;
      iMatDiffusionTimesMu(13,19)=(617*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(14,0)=(dq[0]*dq[1])/480;
      iMatDiffusionTimesMu(14,1)=(-25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(14,2)=0;
      iMatDiffusionTimesMu(14,3)=(-25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(14,4)=(-25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(14,5)=(dq[0]*dq[1])/480;
      iMatDiffusionTimesMu(14,6)=(-25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(14,7)=(dq[0]*dq[1])/480;
      iMatDiffusionTimesMu(14,8)=(89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(14,9)=0;
      iMatDiffusionTimesMu(14,10)=(89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(14,11)=(89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(14,12)=(-239*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(14,13)=(323*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(14,14)=(-11*dq[0]*dq[1])/240;
      iMatDiffusionTimesMu(14,15)=(653*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(14,16)=(323*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(14,17)=(-61*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(14,18)=(653*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(14,19)=(-239*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(15,0)=(89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(15,1)=0;
      iMatDiffusionTimesMu(15,2)=(89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(15,3)=0;
      iMatDiffusionTimesMu(15,4)=0;
      iMatDiffusionTimesMu(15,5)=(89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(15,6)=0;
      iMatDiffusionTimesMu(15,7)=(89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(15,8)=(-89*dq[0]*dq[1])/2160;
      iMatDiffusionTimesMu(15,9)=(-89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(15,10)=(-89*dq[0]*dq[1])/2160;
      iMatDiffusionTimesMu(15,11)=(-89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(15,12)=(287*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(15,13)=(-11*dq[0]*dq[1])/72;
      iMatDiffusionTimesMu(15,14)=(617*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(15,15)=(-11*dq[0]*dq[1])/45;
      iMatDiffusionTimesMu(15,16)=(-11*dq[0]*dq[1])/90;
      iMatDiffusionTimesMu(15,17)=(287*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(15,18)=(-11*dq[0]*dq[1])/72;
      iMatDiffusionTimesMu(15,19)=(617*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(16,0)=(89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(16,1)=0;
      iMatDiffusionTimesMu(16,2)=(89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(16,3)=0;
      iMatDiffusionTimesMu(16,4)=0;
      iMatDiffusionTimesMu(16,5)=(89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(16,6)=0;
      iMatDiffusionTimesMu(16,7)=(89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(16,8)=(-89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(16,9)=(-89*dq[0]*dq[1])/2160;
      iMatDiffusionTimesMu(16,10)=(-89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(16,11)=(-89*dq[0]*dq[1])/2160;
      iMatDiffusionTimesMu(16,12)=(617*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(16,13)=(-11*dq[0]*dq[1])/72;
      iMatDiffusionTimesMu(16,14)=(287*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(16,15)=(-11*dq[0]*dq[1])/90;
      iMatDiffusionTimesMu(16,16)=(-11*dq[0]*dq[1])/45;
      iMatDiffusionTimesMu(16,17)=(617*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(16,18)=(-11*dq[0]*dq[1])/72;
      iMatDiffusionTimesMu(16,19)=(287*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(17,0)=(dq[0]*dq[1])/480;
      iMatDiffusionTimesMu(17,1)=(-25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(17,2)=(dq[0]*dq[1])/480;
      iMatDiffusionTimesMu(17,3)=(-25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(17,4)=(-25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(17,5)=0;
      iMatDiffusionTimesMu(17,6)=(-25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(17,7)=(dq[0]*dq[1])/480;
      iMatDiffusionTimesMu(17,8)=(89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(17,9)=(89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(17,10)=0;
      iMatDiffusionTimesMu(17,11)=(89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(17,12)=(-239*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(17,13)=(653*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(17,14)=(-61*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(17,15)=(323*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(17,16)=(653*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(17,17)=(-11*dq[0]*dq[1])/240;
      iMatDiffusionTimesMu(17,18)=(323*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(17,19)=(-239*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(18,0)=(89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(18,1)=0;
      iMatDiffusionTimesMu(18,2)=(89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(18,3)=0;
      iMatDiffusionTimesMu(18,4)=0;
      iMatDiffusionTimesMu(18,5)=(89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(18,6)=0;
      iMatDiffusionTimesMu(18,7)=(89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(18,8)=(-89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(18,9)=(-89*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(18,10)=(-89*dq[0]*dq[1])/2160;
      iMatDiffusionTimesMu(18,11)=(-89*dq[0]*dq[1])/2160;
      iMatDiffusionTimesMu(18,12)=(617*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(18,13)=(-11*dq[0]*dq[1])/90;
      iMatDiffusionTimesMu(18,14)=(617*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(18,15)=(-11*dq[0]*dq[1])/72;
      iMatDiffusionTimesMu(18,16)=(-11*dq[0]*dq[1])/72;
      iMatDiffusionTimesMu(18,17)=(287*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(18,18)=(-11*dq[0]*dq[1])/45;
      iMatDiffusionTimesMu(18,19)=(287*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(19,0)=(dq[0]*dq[1])/480;
      iMatDiffusionTimesMu(19,1)=(-25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(19,2)=(dq[0]*dq[1])/480;
      iMatDiffusionTimesMu(19,3)=(-25*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(19,4)=(-25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(19,5)=(dq[0]*dq[1])/480;
      iMatDiffusionTimesMu(19,6)=(-25*dq[0]*dq[1])/864;
      iMatDiffusionTimesMu(19,7)=0;
      iMatDiffusionTimesMu(19,8)=(89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(19,9)=(89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(19,10)=(89*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(19,11)=0;
      iMatDiffusionTimesMu(19,12)=(-61*dq[0]*dq[1])/1728;
      iMatDiffusionTimesMu(19,13)=(653*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(19,14)=(-239*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(19,15)=(653*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(19,16)=(323*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(19,17)=(-239*dq[0]*dq[1])/8640;
      iMatDiffusionTimesMu(19,18)=(323*dq[0]*dq[1])/4320;
      iMatDiffusionTimesMu(19,19)=(-11*dq[0]*dq[1])/240;
      lowerMatDiffusionTimesMu(0,0)=(-23*dq[0]*dq[1])/1280;
      lowerMatDiffusionTimesMu(0,1)=(49*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(0,2)=(-5*dq[0]*dq[1])/256;
      lowerMatDiffusionTimesMu(0,3)=(49*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(0,4)=(73*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(0,5)=(-5*dq[0]*dq[1])/256;
      lowerMatDiffusionTimesMu(0,6)=(73*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(0,7)=(-419*dq[0]*dq[1])/23040;
      lowerMatDiffusionTimesMu(0,8)=(dq[0]*dq[1])/48;
      lowerMatDiffusionTimesMu(0,9)=(181*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(0,10)=(181*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(0,11)=(151*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(0,12)=(-73*dq[0]*dq[1])/1280;
      lowerMatDiffusionTimesMu(0,13)=(1117*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(0,14)=(-101*dq[0]*dq[1])/2304;
      lowerMatDiffusionTimesMu(0,15)=(1117*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(0,16)=(503*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(0,17)=(-101*dq[0]*dq[1])/2304;
      lowerMatDiffusionTimesMu(0,18)=(503*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(0,19)=(-1129*dq[0]*dq[1])/23040;
      lowerMatDiffusionTimesMu(1,0)=(571*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(1,1)=(-29*dq[0]*dq[1])/720;
      lowerMatDiffusionTimesMu(1,2)=(571*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(1,3)=(-29*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(1,4)=(-29*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(1,5)=(179*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(1,6)=(-29*dq[0]*dq[1])/1440;
      lowerMatDiffusionTimesMu(1,7)=(179*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(1,8)=(-121*dq[0]*dq[1])/1440;
      lowerMatDiffusionTimesMu(1,9)=(-121*dq[0]*dq[1])/1440;
      lowerMatDiffusionTimesMu(1,10)=(-121*dq[0]*dq[1])/2880;
      lowerMatDiffusionTimesMu(1,11)=(-121*dq[0]*dq[1])/2880;
      lowerMatDiffusionTimesMu(1,12)=(1021*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(1,13)=(-179*dq[0]*dq[1])/720;
      lowerMatDiffusionTimesMu(1,14)=(1021*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(1,15)=(-179*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(1,16)=(-179*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(1,17)=(479*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(1,18)=(-179*dq[0]*dq[1])/1440;
      lowerMatDiffusionTimesMu(1,19)=(479*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(2,0)=(-5*dq[0]*dq[1])/256;
      lowerMatDiffusionTimesMu(2,1)=(49*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(2,2)=(-23*dq[0]*dq[1])/1280;
      lowerMatDiffusionTimesMu(2,3)=(73*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(2,4)=(49*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(2,5)=(-419*dq[0]*dq[1])/23040;
      lowerMatDiffusionTimesMu(2,6)=(73*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(2,7)=(-5*dq[0]*dq[1])/256;
      lowerMatDiffusionTimesMu(2,8)=(181*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(2,9)=(dq[0]*dq[1])/48;
      lowerMatDiffusionTimesMu(2,10)=(151*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(2,11)=(181*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(2,12)=(-101*dq[0]*dq[1])/2304;
      lowerMatDiffusionTimesMu(2,13)=(1117*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(2,14)=(-73*dq[0]*dq[1])/1280;
      lowerMatDiffusionTimesMu(2,15)=(503*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(2,16)=(1117*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(2,17)=(-1129*dq[0]*dq[1])/23040;
      lowerMatDiffusionTimesMu(2,18)=(503*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(2,19)=(-101*dq[0]*dq[1])/2304;
      lowerMatDiffusionTimesMu(3,0)=(571*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(3,1)=(-29*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(3,2)=(179*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(3,3)=(-29*dq[0]*dq[1])/720;
      lowerMatDiffusionTimesMu(3,4)=(-29*dq[0]*dq[1])/1440;
      lowerMatDiffusionTimesMu(3,5)=(571*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(3,6)=(-29*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(3,7)=(179*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(3,8)=(-121*dq[0]*dq[1])/1440;
      lowerMatDiffusionTimesMu(3,9)=(-121*dq[0]*dq[1])/2880;
      lowerMatDiffusionTimesMu(3,10)=(-121*dq[0]*dq[1])/1440;
      lowerMatDiffusionTimesMu(3,11)=(-121*dq[0]*dq[1])/2880;
      lowerMatDiffusionTimesMu(3,12)=(1021*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(3,13)=(-179*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(3,14)=(479*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(3,15)=(-179*dq[0]*dq[1])/720;
      lowerMatDiffusionTimesMu(3,16)=(-179*dq[0]*dq[1])/1440;
      lowerMatDiffusionTimesMu(3,17)=(1021*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(3,18)=(-179*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(3,19)=(479*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(4,0)=(179*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(4,1)=(-29*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(4,2)=(571*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(4,3)=(-29*dq[0]*dq[1])/1440;
      lowerMatDiffusionTimesMu(4,4)=(-29*dq[0]*dq[1])/720;
      lowerMatDiffusionTimesMu(4,5)=(179*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(4,6)=(-29*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(4,7)=(571*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(4,8)=(-121*dq[0]*dq[1])/2880;
      lowerMatDiffusionTimesMu(4,9)=(-121*dq[0]*dq[1])/1440;
      lowerMatDiffusionTimesMu(4,10)=(-121*dq[0]*dq[1])/2880;
      lowerMatDiffusionTimesMu(4,11)=(-121*dq[0]*dq[1])/1440;
      lowerMatDiffusionTimesMu(4,12)=(479*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(4,13)=(-179*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(4,14)=(1021*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(4,15)=(-179*dq[0]*dq[1])/1440;
      lowerMatDiffusionTimesMu(4,16)=(-179*dq[0]*dq[1])/720;
      lowerMatDiffusionTimesMu(4,17)=(479*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(4,18)=(-179*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(4,19)=(1021*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(5,0)=(-5*dq[0]*dq[1])/256;
      lowerMatDiffusionTimesMu(5,1)=(73*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(5,2)=(-419*dq[0]*dq[1])/23040;
      lowerMatDiffusionTimesMu(5,3)=(49*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(5,4)=(73*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(5,5)=(-23*dq[0]*dq[1])/1280;
      lowerMatDiffusionTimesMu(5,6)=(49*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(5,7)=(-5*dq[0]*dq[1])/256;
      lowerMatDiffusionTimesMu(5,8)=(181*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(5,9)=(151*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(5,10)=(dq[0]*dq[1])/48;
      lowerMatDiffusionTimesMu(5,11)=(181*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(5,12)=(-101*dq[0]*dq[1])/2304;
      lowerMatDiffusionTimesMu(5,13)=(503*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(5,14)=(-1129*dq[0]*dq[1])/23040;
      lowerMatDiffusionTimesMu(5,15)=(1117*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(5,16)=(503*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(5,17)=(-73*dq[0]*dq[1])/1280;
      lowerMatDiffusionTimesMu(5,18)=(1117*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(5,19)=(-101*dq[0]*dq[1])/2304;
      lowerMatDiffusionTimesMu(6,0)=(179*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(6,1)=(-29*dq[0]*dq[1])/1440;
      lowerMatDiffusionTimesMu(6,2)=(179*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(6,3)=(-29*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(6,4)=(-29*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(6,5)=(571*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(6,6)=(-29*dq[0]*dq[1])/720;
      lowerMatDiffusionTimesMu(6,7)=(571*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(6,8)=(-121*dq[0]*dq[1])/2880;
      lowerMatDiffusionTimesMu(6,9)=(-121*dq[0]*dq[1])/2880;
      lowerMatDiffusionTimesMu(6,10)=(-121*dq[0]*dq[1])/1440;
      lowerMatDiffusionTimesMu(6,11)=(-121*dq[0]*dq[1])/1440;
      lowerMatDiffusionTimesMu(6,12)=(479*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(6,13)=(-179*dq[0]*dq[1])/1440;
      lowerMatDiffusionTimesMu(6,14)=(479*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(6,15)=(-179*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(6,16)=(-179*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(6,17)=(1021*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(6,18)=(-179*dq[0]*dq[1])/720;
      lowerMatDiffusionTimesMu(6,19)=(1021*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(7,0)=(-419*dq[0]*dq[1])/23040;
      lowerMatDiffusionTimesMu(7,1)=(73*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(7,2)=(-5*dq[0]*dq[1])/256;
      lowerMatDiffusionTimesMu(7,3)=(73*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(7,4)=(49*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(7,5)=(-5*dq[0]*dq[1])/256;
      lowerMatDiffusionTimesMu(7,6)=(49*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(7,7)=(-23*dq[0]*dq[1])/1280;
      lowerMatDiffusionTimesMu(7,8)=(151*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(7,9)=(181*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(7,10)=(181*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(7,11)=(dq[0]*dq[1])/48;
      lowerMatDiffusionTimesMu(7,12)=(-1129*dq[0]*dq[1])/23040;
      lowerMatDiffusionTimesMu(7,13)=(503*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(7,14)=(-101*dq[0]*dq[1])/2304;
      lowerMatDiffusionTimesMu(7,15)=(503*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(7,16)=(1117*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(7,17)=(-101*dq[0]*dq[1])/2304;
      lowerMatDiffusionTimesMu(7,18)=(1117*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(7,19)=(-73*dq[0]*dq[1])/1280;
      lowerMatDiffusionTimesMu(8,0)=(dq[0]*dq[1])/48;
      lowerMatDiffusionTimesMu(8,1)=-(dq[0]*dq[1])/96;
      lowerMatDiffusionTimesMu(8,2)=(5*dq[0]*dq[1])/384;
      lowerMatDiffusionTimesMu(8,3)=-(dq[0]*dq[1])/96;
      lowerMatDiffusionTimesMu(8,4)=-(dq[0]*dq[1])/192;
      lowerMatDiffusionTimesMu(8,5)=(5*dq[0]*dq[1])/384;
      lowerMatDiffusionTimesMu(8,6)=-(dq[0]*dq[1])/192;
      lowerMatDiffusionTimesMu(8,7)=(dq[0]*dq[1])/128;
      lowerMatDiffusionTimesMu(8,8)=-(dq[0]*dq[1])/24;
      lowerMatDiffusionTimesMu(8,9)=-(dq[0]*dq[1])/48;
      lowerMatDiffusionTimesMu(8,10)=-(dq[0]*dq[1])/48;
      lowerMatDiffusionTimesMu(8,11)=-(dq[0]*dq[1])/96;
      lowerMatDiffusionTimesMu(8,12)=(dq[0]*dq[1])/48;
      lowerMatDiffusionTimesMu(8,13)=(-29*dq[0]*dq[1])/288;
      lowerMatDiffusionTimesMu(8,14)=(41*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(8,15)=(-29*dq[0]*dq[1])/288;
      lowerMatDiffusionTimesMu(8,16)=(-29*dq[0]*dq[1])/576;
      lowerMatDiffusionTimesMu(8,17)=(41*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(8,18)=(-29*dq[0]*dq[1])/576;
      lowerMatDiffusionTimesMu(8,19)=(35*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(9,0)=(5*dq[0]*dq[1])/384;
      lowerMatDiffusionTimesMu(9,1)=-(dq[0]*dq[1])/96;
      lowerMatDiffusionTimesMu(9,2)=(dq[0]*dq[1])/48;
      lowerMatDiffusionTimesMu(9,3)=-(dq[0]*dq[1])/192;
      lowerMatDiffusionTimesMu(9,4)=-(dq[0]*dq[1])/96;
      lowerMatDiffusionTimesMu(9,5)=(dq[0]*dq[1])/128;
      lowerMatDiffusionTimesMu(9,6)=-(dq[0]*dq[1])/192;
      lowerMatDiffusionTimesMu(9,7)=(5*dq[0]*dq[1])/384;
      lowerMatDiffusionTimesMu(9,8)=-(dq[0]*dq[1])/48;
      lowerMatDiffusionTimesMu(9,9)=-(dq[0]*dq[1])/24;
      lowerMatDiffusionTimesMu(9,10)=-(dq[0]*dq[1])/96;
      lowerMatDiffusionTimesMu(9,11)=-(dq[0]*dq[1])/48;
      lowerMatDiffusionTimesMu(9,12)=(41*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(9,13)=(-29*dq[0]*dq[1])/288;
      lowerMatDiffusionTimesMu(9,14)=(dq[0]*dq[1])/48;
      lowerMatDiffusionTimesMu(9,15)=(-29*dq[0]*dq[1])/576;
      lowerMatDiffusionTimesMu(9,16)=(-29*dq[0]*dq[1])/288;
      lowerMatDiffusionTimesMu(9,17)=(35*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(9,18)=(-29*dq[0]*dq[1])/576;
      lowerMatDiffusionTimesMu(9,19)=(41*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(10,0)=(5*dq[0]*dq[1])/384;
      lowerMatDiffusionTimesMu(10,1)=-(dq[0]*dq[1])/192;
      lowerMatDiffusionTimesMu(10,2)=(dq[0]*dq[1])/128;
      lowerMatDiffusionTimesMu(10,3)=-(dq[0]*dq[1])/96;
      lowerMatDiffusionTimesMu(10,4)=-(dq[0]*dq[1])/192;
      lowerMatDiffusionTimesMu(10,5)=(dq[0]*dq[1])/48;
      lowerMatDiffusionTimesMu(10,6)=-(dq[0]*dq[1])/96;
      lowerMatDiffusionTimesMu(10,7)=(5*dq[0]*dq[1])/384;
      lowerMatDiffusionTimesMu(10,8)=-(dq[0]*dq[1])/48;
      lowerMatDiffusionTimesMu(10,9)=-(dq[0]*dq[1])/96;
      lowerMatDiffusionTimesMu(10,10)=-(dq[0]*dq[1])/24;
      lowerMatDiffusionTimesMu(10,11)=-(dq[0]*dq[1])/48;
      lowerMatDiffusionTimesMu(10,12)=(41*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(10,13)=(-29*dq[0]*dq[1])/576;
      lowerMatDiffusionTimesMu(10,14)=(35*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(10,15)=(-29*dq[0]*dq[1])/288;
      lowerMatDiffusionTimesMu(10,16)=(-29*dq[0]*dq[1])/576;
      lowerMatDiffusionTimesMu(10,17)=(dq[0]*dq[1])/48;
      lowerMatDiffusionTimesMu(10,18)=(-29*dq[0]*dq[1])/288;
      lowerMatDiffusionTimesMu(10,19)=(41*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(11,0)=(dq[0]*dq[1])/128;
      lowerMatDiffusionTimesMu(11,1)=-(dq[0]*dq[1])/192;
      lowerMatDiffusionTimesMu(11,2)=(5*dq[0]*dq[1])/384;
      lowerMatDiffusionTimesMu(11,3)=-(dq[0]*dq[1])/192;
      lowerMatDiffusionTimesMu(11,4)=-(dq[0]*dq[1])/96;
      lowerMatDiffusionTimesMu(11,5)=(5*dq[0]*dq[1])/384;
      lowerMatDiffusionTimesMu(11,6)=-(dq[0]*dq[1])/96;
      lowerMatDiffusionTimesMu(11,7)=(dq[0]*dq[1])/48;
      lowerMatDiffusionTimesMu(11,8)=-(dq[0]*dq[1])/96;
      lowerMatDiffusionTimesMu(11,9)=-(dq[0]*dq[1])/48;
      lowerMatDiffusionTimesMu(11,10)=-(dq[0]*dq[1])/48;
      lowerMatDiffusionTimesMu(11,11)=-(dq[0]*dq[1])/24;
      lowerMatDiffusionTimesMu(11,12)=(35*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(11,13)=(-29*dq[0]*dq[1])/576;
      lowerMatDiffusionTimesMu(11,14)=(41*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(11,15)=(-29*dq[0]*dq[1])/576;
      lowerMatDiffusionTimesMu(11,16)=(-29*dq[0]*dq[1])/288;
      lowerMatDiffusionTimesMu(11,17)=(41*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(11,18)=(-29*dq[0]*dq[1])/288;
      lowerMatDiffusionTimesMu(11,19)=(dq[0]*dq[1])/48;
      lowerMatDiffusionTimesMu(12,0)=(-43*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(12,1)=(23*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(12,2)=(-31*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(12,3)=(23*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(12,4)=(7*dq[0]*dq[1])/1920;
      lowerMatDiffusionTimesMu(12,5)=(-31*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(12,6)=(7*dq[0]*dq[1])/1920;
      lowerMatDiffusionTimesMu(12,7)=(-43*dq[0]*dq[1])/7680;
      lowerMatDiffusionTimesMu(12,8)=(dq[0]*dq[1])/48;
      lowerMatDiffusionTimesMu(12,9)=(5*dq[0]*dq[1])/384;
      lowerMatDiffusionTimesMu(12,10)=(5*dq[0]*dq[1])/384;
      lowerMatDiffusionTimesMu(12,11)=(dq[0]*dq[1])/128;
      lowerMatDiffusionTimesMu(12,12)=(-23*dq[0]*dq[1])/1280;
      lowerMatDiffusionTimesMu(12,13)=(667*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(12,14)=(-83*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(12,15)=(667*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(12,16)=(203*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(12,17)=(-83*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(12,18)=(203*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(12,19)=(-467*dq[0]*dq[1])/23040;
      lowerMatDiffusionTimesMu(13,0)=(23*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(13,1)=-(dq[0]*dq[1])/240;
      lowerMatDiffusionTimesMu(13,2)=(23*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(13,3)=-(dq[0]*dq[1])/384;
      lowerMatDiffusionTimesMu(13,4)=-(dq[0]*dq[1])/384;
      lowerMatDiffusionTimesMu(13,5)=(7*dq[0]*dq[1])/1920;
      lowerMatDiffusionTimesMu(13,6)=-(dq[0]*dq[1])/480;
      lowerMatDiffusionTimesMu(13,7)=(7*dq[0]*dq[1])/1920;
      lowerMatDiffusionTimesMu(13,8)=-(dq[0]*dq[1])/96;
      lowerMatDiffusionTimesMu(13,9)=-(dq[0]*dq[1])/96;
      lowerMatDiffusionTimesMu(13,10)=-(dq[0]*dq[1])/192;
      lowerMatDiffusionTimesMu(13,11)=-(dq[0]*dq[1])/192;
      lowerMatDiffusionTimesMu(13,12)=(49*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(13,13)=(-29*dq[0]*dq[1])/720;
      lowerMatDiffusionTimesMu(13,14)=(49*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(13,15)=(-29*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(13,16)=(-29*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(13,17)=(73*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(13,18)=(-29*dq[0]*dq[1])/1440;
      lowerMatDiffusionTimesMu(13,19)=(73*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(14,0)=(-31*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(14,1)=(23*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(14,2)=(-43*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(14,3)=(7*dq[0]*dq[1])/1920;
      lowerMatDiffusionTimesMu(14,4)=(23*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(14,5)=(-43*dq[0]*dq[1])/7680;
      lowerMatDiffusionTimesMu(14,6)=(7*dq[0]*dq[1])/1920;
      lowerMatDiffusionTimesMu(14,7)=(-31*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(14,8)=(5*dq[0]*dq[1])/384;
      lowerMatDiffusionTimesMu(14,9)=(dq[0]*dq[1])/48;
      lowerMatDiffusionTimesMu(14,10)=(dq[0]*dq[1])/128;
      lowerMatDiffusionTimesMu(14,11)=(5*dq[0]*dq[1])/384;
      lowerMatDiffusionTimesMu(14,12)=(-83*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(14,13)=(667*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(14,14)=(-23*dq[0]*dq[1])/1280;
      lowerMatDiffusionTimesMu(14,15)=(203*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(14,16)=(667*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(14,17)=(-467*dq[0]*dq[1])/23040;
      lowerMatDiffusionTimesMu(14,18)=(203*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(14,19)=(-83*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(15,0)=(23*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(15,1)=-(dq[0]*dq[1])/384;
      lowerMatDiffusionTimesMu(15,2)=(7*dq[0]*dq[1])/1920;
      lowerMatDiffusionTimesMu(15,3)=-(dq[0]*dq[1])/240;
      lowerMatDiffusionTimesMu(15,4)=-(dq[0]*dq[1])/480;
      lowerMatDiffusionTimesMu(15,5)=(23*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(15,6)=-(dq[0]*dq[1])/384;
      lowerMatDiffusionTimesMu(15,7)=(7*dq[0]*dq[1])/1920;
      lowerMatDiffusionTimesMu(15,8)=-(dq[0]*dq[1])/96;
      lowerMatDiffusionTimesMu(15,9)=-(dq[0]*dq[1])/192;
      lowerMatDiffusionTimesMu(15,10)=-(dq[0]*dq[1])/96;
      lowerMatDiffusionTimesMu(15,11)=-(dq[0]*dq[1])/192;
      lowerMatDiffusionTimesMu(15,12)=(49*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(15,13)=(-29*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(15,14)=(73*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(15,15)=(-29*dq[0]*dq[1])/720;
      lowerMatDiffusionTimesMu(15,16)=(-29*dq[0]*dq[1])/1440;
      lowerMatDiffusionTimesMu(15,17)=(49*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(15,18)=(-29*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(15,19)=(73*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(16,0)=(7*dq[0]*dq[1])/1920;
      lowerMatDiffusionTimesMu(16,1)=-(dq[0]*dq[1])/384;
      lowerMatDiffusionTimesMu(16,2)=(23*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(16,3)=-(dq[0]*dq[1])/480;
      lowerMatDiffusionTimesMu(16,4)=-(dq[0]*dq[1])/240;
      lowerMatDiffusionTimesMu(16,5)=(7*dq[0]*dq[1])/1920;
      lowerMatDiffusionTimesMu(16,6)=-(dq[0]*dq[1])/384;
      lowerMatDiffusionTimesMu(16,7)=(23*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(16,8)=-(dq[0]*dq[1])/192;
      lowerMatDiffusionTimesMu(16,9)=-(dq[0]*dq[1])/96;
      lowerMatDiffusionTimesMu(16,10)=-(dq[0]*dq[1])/192;
      lowerMatDiffusionTimesMu(16,11)=-(dq[0]*dq[1])/96;
      lowerMatDiffusionTimesMu(16,12)=(73*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(16,13)=(-29*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(16,14)=(49*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(16,15)=(-29*dq[0]*dq[1])/1440;
      lowerMatDiffusionTimesMu(16,16)=(-29*dq[0]*dq[1])/720;
      lowerMatDiffusionTimesMu(16,17)=(73*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(16,18)=(-29*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(16,19)=(49*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(17,0)=(-31*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(17,1)=(7*dq[0]*dq[1])/1920;
      lowerMatDiffusionTimesMu(17,2)=(-43*dq[0]*dq[1])/7680;
      lowerMatDiffusionTimesMu(17,3)=(23*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(17,4)=(7*dq[0]*dq[1])/1920;
      lowerMatDiffusionTimesMu(17,5)=(-43*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(17,6)=(23*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(17,7)=(-31*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(17,8)=(5*dq[0]*dq[1])/384;
      lowerMatDiffusionTimesMu(17,9)=(dq[0]*dq[1])/128;
      lowerMatDiffusionTimesMu(17,10)=(dq[0]*dq[1])/48;
      lowerMatDiffusionTimesMu(17,11)=(5*dq[0]*dq[1])/384;
      lowerMatDiffusionTimesMu(17,12)=(-83*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(17,13)=(203*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(17,14)=(-467*dq[0]*dq[1])/23040;
      lowerMatDiffusionTimesMu(17,15)=(667*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(17,16)=(203*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(17,17)=(-23*dq[0]*dq[1])/1280;
      lowerMatDiffusionTimesMu(17,18)=(667*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(17,19)=(-83*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(18,0)=(7*dq[0]*dq[1])/1920;
      lowerMatDiffusionTimesMu(18,1)=-(dq[0]*dq[1])/480;
      lowerMatDiffusionTimesMu(18,2)=(7*dq[0]*dq[1])/1920;
      lowerMatDiffusionTimesMu(18,3)=-(dq[0]*dq[1])/384;
      lowerMatDiffusionTimesMu(18,4)=-(dq[0]*dq[1])/384;
      lowerMatDiffusionTimesMu(18,5)=(23*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(18,6)=-(dq[0]*dq[1])/240;
      lowerMatDiffusionTimesMu(18,7)=(23*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(18,8)=-(dq[0]*dq[1])/192;
      lowerMatDiffusionTimesMu(18,9)=-(dq[0]*dq[1])/192;
      lowerMatDiffusionTimesMu(18,10)=-(dq[0]*dq[1])/96;
      lowerMatDiffusionTimesMu(18,11)=-(dq[0]*dq[1])/96;
      lowerMatDiffusionTimesMu(18,12)=(73*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(18,13)=(-29*dq[0]*dq[1])/1440;
      lowerMatDiffusionTimesMu(18,14)=(73*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(18,15)=(-29*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(18,16)=(-29*dq[0]*dq[1])/1152;
      lowerMatDiffusionTimesMu(18,17)=(49*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(18,18)=(-29*dq[0]*dq[1])/720;
      lowerMatDiffusionTimesMu(18,19)=(49*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(19,0)=(-43*dq[0]*dq[1])/7680;
      lowerMatDiffusionTimesMu(19,1)=(7*dq[0]*dq[1])/1920;
      lowerMatDiffusionTimesMu(19,2)=(-31*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(19,3)=(7*dq[0]*dq[1])/1920;
      lowerMatDiffusionTimesMu(19,4)=(23*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(19,5)=(-31*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(19,6)=(23*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(19,7)=(-43*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(19,8)=(dq[0]*dq[1])/128;
      lowerMatDiffusionTimesMu(19,9)=(5*dq[0]*dq[1])/384;
      lowerMatDiffusionTimesMu(19,10)=(5*dq[0]*dq[1])/384;
      lowerMatDiffusionTimesMu(19,11)=(dq[0]*dq[1])/48;
      lowerMatDiffusionTimesMu(19,12)=(-467*dq[0]*dq[1])/23040;
      lowerMatDiffusionTimesMu(19,13)=(203*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(19,14)=(-83*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(19,15)=(203*dq[0]*dq[1])/5760;
      lowerMatDiffusionTimesMu(19,16)=(667*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(19,17)=(-83*dq[0]*dq[1])/3840;
      lowerMatDiffusionTimesMu(19,18)=(667*dq[0]*dq[1])/11520;
      lowerMatDiffusionTimesMu(19,19)=(-23*dq[0]*dq[1])/1280;
      upperMatDiffusionTimesMu(0,0)=(23*dq[0]*dq[1])/1280;
      upperMatDiffusionTimesMu(0,1)=(-667*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(0,2)=(83*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(0,3)=(-667*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(0,4)=(-203*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(0,5)=(83*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(0,6)=(-203*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(0,7)=(467*dq[0]*dq[1])/23040;
      upperMatDiffusionTimesMu(0,8)=-(dq[0]*dq[1])/48;
      upperMatDiffusionTimesMu(0,9)=(-5*dq[0]*dq[1])/384;
      upperMatDiffusionTimesMu(0,10)=(-5*dq[0]*dq[1])/384;
      upperMatDiffusionTimesMu(0,11)=-(dq[0]*dq[1])/128;
      upperMatDiffusionTimesMu(0,12)=(43*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(0,13)=(-23*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(0,14)=(31*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(0,15)=(-23*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(0,16)=(-7*dq[0]*dq[1])/1920;
      upperMatDiffusionTimesMu(0,17)=(31*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(0,18)=(-7*dq[0]*dq[1])/1920;
      upperMatDiffusionTimesMu(0,19)=(43*dq[0]*dq[1])/7680;
      upperMatDiffusionTimesMu(1,0)=(-49*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(1,1)=(29*dq[0]*dq[1])/720;
      upperMatDiffusionTimesMu(1,2)=(-49*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(1,3)=(29*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(1,4)=(29*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(1,5)=(-73*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(1,6)=(29*dq[0]*dq[1])/1440;
      upperMatDiffusionTimesMu(1,7)=(-73*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(1,8)=(dq[0]*dq[1])/96;
      upperMatDiffusionTimesMu(1,9)=(dq[0]*dq[1])/96;
      upperMatDiffusionTimesMu(1,10)=(dq[0]*dq[1])/192;
      upperMatDiffusionTimesMu(1,11)=(dq[0]*dq[1])/192;
      upperMatDiffusionTimesMu(1,12)=(-23*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(1,13)=(dq[0]*dq[1])/240;
      upperMatDiffusionTimesMu(1,14)=(-23*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(1,15)=(dq[0]*dq[1])/384;
      upperMatDiffusionTimesMu(1,16)=(dq[0]*dq[1])/384;
      upperMatDiffusionTimesMu(1,17)=(-7*dq[0]*dq[1])/1920;
      upperMatDiffusionTimesMu(1,18)=(dq[0]*dq[1])/480;
      upperMatDiffusionTimesMu(1,19)=(-7*dq[0]*dq[1])/1920;
      upperMatDiffusionTimesMu(2,0)=(83*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(2,1)=(-667*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(2,2)=(23*dq[0]*dq[1])/1280;
      upperMatDiffusionTimesMu(2,3)=(-203*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(2,4)=(-667*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(2,5)=(467*dq[0]*dq[1])/23040;
      upperMatDiffusionTimesMu(2,6)=(-203*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(2,7)=(83*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(2,8)=(-5*dq[0]*dq[1])/384;
      upperMatDiffusionTimesMu(2,9)=-(dq[0]*dq[1])/48;
      upperMatDiffusionTimesMu(2,10)=-(dq[0]*dq[1])/128;
      upperMatDiffusionTimesMu(2,11)=(-5*dq[0]*dq[1])/384;
      upperMatDiffusionTimesMu(2,12)=(31*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(2,13)=(-23*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(2,14)=(43*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(2,15)=(-7*dq[0]*dq[1])/1920;
      upperMatDiffusionTimesMu(2,16)=(-23*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(2,17)=(43*dq[0]*dq[1])/7680;
      upperMatDiffusionTimesMu(2,18)=(-7*dq[0]*dq[1])/1920;
      upperMatDiffusionTimesMu(2,19)=(31*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(3,0)=(-49*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(3,1)=(29*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(3,2)=(-73*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(3,3)=(29*dq[0]*dq[1])/720;
      upperMatDiffusionTimesMu(3,4)=(29*dq[0]*dq[1])/1440;
      upperMatDiffusionTimesMu(3,5)=(-49*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(3,6)=(29*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(3,7)=(-73*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(3,8)=(dq[0]*dq[1])/96;
      upperMatDiffusionTimesMu(3,9)=(dq[0]*dq[1])/192;
      upperMatDiffusionTimesMu(3,10)=(dq[0]*dq[1])/96;
      upperMatDiffusionTimesMu(3,11)=(dq[0]*dq[1])/192;
      upperMatDiffusionTimesMu(3,12)=(-23*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(3,13)=(dq[0]*dq[1])/384;
      upperMatDiffusionTimesMu(3,14)=(-7*dq[0]*dq[1])/1920;
      upperMatDiffusionTimesMu(3,15)=(dq[0]*dq[1])/240;
      upperMatDiffusionTimesMu(3,16)=(dq[0]*dq[1])/480;
      upperMatDiffusionTimesMu(3,17)=(-23*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(3,18)=(dq[0]*dq[1])/384;
      upperMatDiffusionTimesMu(3,19)=(-7*dq[0]*dq[1])/1920;
      upperMatDiffusionTimesMu(4,0)=(-73*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(4,1)=(29*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(4,2)=(-49*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(4,3)=(29*dq[0]*dq[1])/1440;
      upperMatDiffusionTimesMu(4,4)=(29*dq[0]*dq[1])/720;
      upperMatDiffusionTimesMu(4,5)=(-73*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(4,6)=(29*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(4,7)=(-49*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(4,8)=(dq[0]*dq[1])/192;
      upperMatDiffusionTimesMu(4,9)=(dq[0]*dq[1])/96;
      upperMatDiffusionTimesMu(4,10)=(dq[0]*dq[1])/192;
      upperMatDiffusionTimesMu(4,11)=(dq[0]*dq[1])/96;
      upperMatDiffusionTimesMu(4,12)=(-7*dq[0]*dq[1])/1920;
      upperMatDiffusionTimesMu(4,13)=(dq[0]*dq[1])/384;
      upperMatDiffusionTimesMu(4,14)=(-23*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(4,15)=(dq[0]*dq[1])/480;
      upperMatDiffusionTimesMu(4,16)=(dq[0]*dq[1])/240;
      upperMatDiffusionTimesMu(4,17)=(-7*dq[0]*dq[1])/1920;
      upperMatDiffusionTimesMu(4,18)=(dq[0]*dq[1])/384;
      upperMatDiffusionTimesMu(4,19)=(-23*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(5,0)=(83*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(5,1)=(-203*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(5,2)=(467*dq[0]*dq[1])/23040;
      upperMatDiffusionTimesMu(5,3)=(-667*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(5,4)=(-203*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(5,5)=(23*dq[0]*dq[1])/1280;
      upperMatDiffusionTimesMu(5,6)=(-667*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(5,7)=(83*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(5,8)=(-5*dq[0]*dq[1])/384;
      upperMatDiffusionTimesMu(5,9)=-(dq[0]*dq[1])/128;
      upperMatDiffusionTimesMu(5,10)=-(dq[0]*dq[1])/48;
      upperMatDiffusionTimesMu(5,11)=(-5*dq[0]*dq[1])/384;
      upperMatDiffusionTimesMu(5,12)=(31*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(5,13)=(-7*dq[0]*dq[1])/1920;
      upperMatDiffusionTimesMu(5,14)=(43*dq[0]*dq[1])/7680;
      upperMatDiffusionTimesMu(5,15)=(-23*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(5,16)=(-7*dq[0]*dq[1])/1920;
      upperMatDiffusionTimesMu(5,17)=(43*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(5,18)=(-23*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(5,19)=(31*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(6,0)=(-73*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(6,1)=(29*dq[0]*dq[1])/1440;
      upperMatDiffusionTimesMu(6,2)=(-73*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(6,3)=(29*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(6,4)=(29*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(6,5)=(-49*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(6,6)=(29*dq[0]*dq[1])/720;
      upperMatDiffusionTimesMu(6,7)=(-49*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(6,8)=(dq[0]*dq[1])/192;
      upperMatDiffusionTimesMu(6,9)=(dq[0]*dq[1])/192;
      upperMatDiffusionTimesMu(6,10)=(dq[0]*dq[1])/96;
      upperMatDiffusionTimesMu(6,11)=(dq[0]*dq[1])/96;
      upperMatDiffusionTimesMu(6,12)=(-7*dq[0]*dq[1])/1920;
      upperMatDiffusionTimesMu(6,13)=(dq[0]*dq[1])/480;
      upperMatDiffusionTimesMu(6,14)=(-7*dq[0]*dq[1])/1920;
      upperMatDiffusionTimesMu(6,15)=(dq[0]*dq[1])/384;
      upperMatDiffusionTimesMu(6,16)=(dq[0]*dq[1])/384;
      upperMatDiffusionTimesMu(6,17)=(-23*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(6,18)=(dq[0]*dq[1])/240;
      upperMatDiffusionTimesMu(6,19)=(-23*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(7,0)=(467*dq[0]*dq[1])/23040;
      upperMatDiffusionTimesMu(7,1)=(-203*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(7,2)=(83*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(7,3)=(-203*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(7,4)=(-667*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(7,5)=(83*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(7,6)=(-667*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(7,7)=(23*dq[0]*dq[1])/1280;
      upperMatDiffusionTimesMu(7,8)=-(dq[0]*dq[1])/128;
      upperMatDiffusionTimesMu(7,9)=(-5*dq[0]*dq[1])/384;
      upperMatDiffusionTimesMu(7,10)=(-5*dq[0]*dq[1])/384;
      upperMatDiffusionTimesMu(7,11)=-(dq[0]*dq[1])/48;
      upperMatDiffusionTimesMu(7,12)=(43*dq[0]*dq[1])/7680;
      upperMatDiffusionTimesMu(7,13)=(-7*dq[0]*dq[1])/1920;
      upperMatDiffusionTimesMu(7,14)=(31*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(7,15)=(-7*dq[0]*dq[1])/1920;
      upperMatDiffusionTimesMu(7,16)=(-23*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(7,17)=(31*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(7,18)=(-23*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(7,19)=(43*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(8,0)=-(dq[0]*dq[1])/48;
      upperMatDiffusionTimesMu(8,1)=(29*dq[0]*dq[1])/288;
      upperMatDiffusionTimesMu(8,2)=(-41*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(8,3)=(29*dq[0]*dq[1])/288;
      upperMatDiffusionTimesMu(8,4)=(29*dq[0]*dq[1])/576;
      upperMatDiffusionTimesMu(8,5)=(-41*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(8,6)=(29*dq[0]*dq[1])/576;
      upperMatDiffusionTimesMu(8,7)=(-35*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(8,8)=(dq[0]*dq[1])/24;
      upperMatDiffusionTimesMu(8,9)=(dq[0]*dq[1])/48;
      upperMatDiffusionTimesMu(8,10)=(dq[0]*dq[1])/48;
      upperMatDiffusionTimesMu(8,11)=(dq[0]*dq[1])/96;
      upperMatDiffusionTimesMu(8,12)=-(dq[0]*dq[1])/48;
      upperMatDiffusionTimesMu(8,13)=(dq[0]*dq[1])/96;
      upperMatDiffusionTimesMu(8,14)=(-5*dq[0]*dq[1])/384;
      upperMatDiffusionTimesMu(8,15)=(dq[0]*dq[1])/96;
      upperMatDiffusionTimesMu(8,16)=(dq[0]*dq[1])/192;
      upperMatDiffusionTimesMu(8,17)=(-5*dq[0]*dq[1])/384;
      upperMatDiffusionTimesMu(8,18)=(dq[0]*dq[1])/192;
      upperMatDiffusionTimesMu(8,19)=-(dq[0]*dq[1])/128;
      upperMatDiffusionTimesMu(9,0)=(-41*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(9,1)=(29*dq[0]*dq[1])/288;
      upperMatDiffusionTimesMu(9,2)=-(dq[0]*dq[1])/48;
      upperMatDiffusionTimesMu(9,3)=(29*dq[0]*dq[1])/576;
      upperMatDiffusionTimesMu(9,4)=(29*dq[0]*dq[1])/288;
      upperMatDiffusionTimesMu(9,5)=(-35*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(9,6)=(29*dq[0]*dq[1])/576;
      upperMatDiffusionTimesMu(9,7)=(-41*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(9,8)=(dq[0]*dq[1])/48;
      upperMatDiffusionTimesMu(9,9)=(dq[0]*dq[1])/24;
      upperMatDiffusionTimesMu(9,10)=(dq[0]*dq[1])/96;
      upperMatDiffusionTimesMu(9,11)=(dq[0]*dq[1])/48;
      upperMatDiffusionTimesMu(9,12)=(-5*dq[0]*dq[1])/384;
      upperMatDiffusionTimesMu(9,13)=(dq[0]*dq[1])/96;
      upperMatDiffusionTimesMu(9,14)=-(dq[0]*dq[1])/48;
      upperMatDiffusionTimesMu(9,15)=(dq[0]*dq[1])/192;
      upperMatDiffusionTimesMu(9,16)=(dq[0]*dq[1])/96;
      upperMatDiffusionTimesMu(9,17)=-(dq[0]*dq[1])/128;
      upperMatDiffusionTimesMu(9,18)=(dq[0]*dq[1])/192;
      upperMatDiffusionTimesMu(9,19)=(-5*dq[0]*dq[1])/384;
      upperMatDiffusionTimesMu(10,0)=(-41*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(10,1)=(29*dq[0]*dq[1])/576;
      upperMatDiffusionTimesMu(10,2)=(-35*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(10,3)=(29*dq[0]*dq[1])/288;
      upperMatDiffusionTimesMu(10,4)=(29*dq[0]*dq[1])/576;
      upperMatDiffusionTimesMu(10,5)=-(dq[0]*dq[1])/48;
      upperMatDiffusionTimesMu(10,6)=(29*dq[0]*dq[1])/288;
      upperMatDiffusionTimesMu(10,7)=(-41*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(10,8)=(dq[0]*dq[1])/48;
      upperMatDiffusionTimesMu(10,9)=(dq[0]*dq[1])/96;
      upperMatDiffusionTimesMu(10,10)=(dq[0]*dq[1])/24;
      upperMatDiffusionTimesMu(10,11)=(dq[0]*dq[1])/48;
      upperMatDiffusionTimesMu(10,12)=(-5*dq[0]*dq[1])/384;
      upperMatDiffusionTimesMu(10,13)=(dq[0]*dq[1])/192;
      upperMatDiffusionTimesMu(10,14)=-(dq[0]*dq[1])/128;
      upperMatDiffusionTimesMu(10,15)=(dq[0]*dq[1])/96;
      upperMatDiffusionTimesMu(10,16)=(dq[0]*dq[1])/192;
      upperMatDiffusionTimesMu(10,17)=-(dq[0]*dq[1])/48;
      upperMatDiffusionTimesMu(10,18)=(dq[0]*dq[1])/96;
      upperMatDiffusionTimesMu(10,19)=(-5*dq[0]*dq[1])/384;
      upperMatDiffusionTimesMu(11,0)=(-35*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(11,1)=(29*dq[0]*dq[1])/576;
      upperMatDiffusionTimesMu(11,2)=(-41*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(11,3)=(29*dq[0]*dq[1])/576;
      upperMatDiffusionTimesMu(11,4)=(29*dq[0]*dq[1])/288;
      upperMatDiffusionTimesMu(11,5)=(-41*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(11,6)=(29*dq[0]*dq[1])/288;
      upperMatDiffusionTimesMu(11,7)=-(dq[0]*dq[1])/48;
      upperMatDiffusionTimesMu(11,8)=(dq[0]*dq[1])/96;
      upperMatDiffusionTimesMu(11,9)=(dq[0]*dq[1])/48;
      upperMatDiffusionTimesMu(11,10)=(dq[0]*dq[1])/48;
      upperMatDiffusionTimesMu(11,11)=(dq[0]*dq[1])/24;
      upperMatDiffusionTimesMu(11,12)=-(dq[0]*dq[1])/128;
      upperMatDiffusionTimesMu(11,13)=(dq[0]*dq[1])/192;
      upperMatDiffusionTimesMu(11,14)=(-5*dq[0]*dq[1])/384;
      upperMatDiffusionTimesMu(11,15)=(dq[0]*dq[1])/192;
      upperMatDiffusionTimesMu(11,16)=(dq[0]*dq[1])/96;
      upperMatDiffusionTimesMu(11,17)=(-5*dq[0]*dq[1])/384;
      upperMatDiffusionTimesMu(11,18)=(dq[0]*dq[1])/96;
      upperMatDiffusionTimesMu(11,19)=-(dq[0]*dq[1])/48;
      upperMatDiffusionTimesMu(12,0)=(73*dq[0]*dq[1])/1280;
      upperMatDiffusionTimesMu(12,1)=(-1117*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(12,2)=(101*dq[0]*dq[1])/2304;
      upperMatDiffusionTimesMu(12,3)=(-1117*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(12,4)=(-503*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(12,5)=(101*dq[0]*dq[1])/2304;
      upperMatDiffusionTimesMu(12,6)=(-503*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(12,7)=(1129*dq[0]*dq[1])/23040;
      upperMatDiffusionTimesMu(12,8)=-(dq[0]*dq[1])/48;
      upperMatDiffusionTimesMu(12,9)=(-181*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(12,10)=(-181*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(12,11)=(-151*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(12,12)=(23*dq[0]*dq[1])/1280;
      upperMatDiffusionTimesMu(12,13)=(-49*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(12,14)=(5*dq[0]*dq[1])/256;
      upperMatDiffusionTimesMu(12,15)=(-49*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(12,16)=(-73*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(12,17)=(5*dq[0]*dq[1])/256;
      upperMatDiffusionTimesMu(12,18)=(-73*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(12,19)=(419*dq[0]*dq[1])/23040;
      upperMatDiffusionTimesMu(13,0)=(-1021*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(13,1)=(179*dq[0]*dq[1])/720;
      upperMatDiffusionTimesMu(13,2)=(-1021*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(13,3)=(179*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(13,4)=(179*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(13,5)=(-479*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(13,6)=(179*dq[0]*dq[1])/1440;
      upperMatDiffusionTimesMu(13,7)=(-479*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(13,8)=(121*dq[0]*dq[1])/1440;
      upperMatDiffusionTimesMu(13,9)=(121*dq[0]*dq[1])/1440;
      upperMatDiffusionTimesMu(13,10)=(121*dq[0]*dq[1])/2880;
      upperMatDiffusionTimesMu(13,11)=(121*dq[0]*dq[1])/2880;
      upperMatDiffusionTimesMu(13,12)=(-571*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(13,13)=(29*dq[0]*dq[1])/720;
      upperMatDiffusionTimesMu(13,14)=(-571*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(13,15)=(29*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(13,16)=(29*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(13,17)=(-179*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(13,18)=(29*dq[0]*dq[1])/1440;
      upperMatDiffusionTimesMu(13,19)=(-179*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(14,0)=(101*dq[0]*dq[1])/2304;
      upperMatDiffusionTimesMu(14,1)=(-1117*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(14,2)=(73*dq[0]*dq[1])/1280;
      upperMatDiffusionTimesMu(14,3)=(-503*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(14,4)=(-1117*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(14,5)=(1129*dq[0]*dq[1])/23040;
      upperMatDiffusionTimesMu(14,6)=(-503*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(14,7)=(101*dq[0]*dq[1])/2304;
      upperMatDiffusionTimesMu(14,8)=(-181*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(14,9)=-(dq[0]*dq[1])/48;
      upperMatDiffusionTimesMu(14,10)=(-151*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(14,11)=(-181*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(14,12)=(5*dq[0]*dq[1])/256;
      upperMatDiffusionTimesMu(14,13)=(-49*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(14,14)=(23*dq[0]*dq[1])/1280;
      upperMatDiffusionTimesMu(14,15)=(-73*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(14,16)=(-49*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(14,17)=(419*dq[0]*dq[1])/23040;
      upperMatDiffusionTimesMu(14,18)=(-73*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(14,19)=(5*dq[0]*dq[1])/256;
      upperMatDiffusionTimesMu(15,0)=(-1021*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(15,1)=(179*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(15,2)=(-479*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(15,3)=(179*dq[0]*dq[1])/720;
      upperMatDiffusionTimesMu(15,4)=(179*dq[0]*dq[1])/1440;
      upperMatDiffusionTimesMu(15,5)=(-1021*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(15,6)=(179*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(15,7)=(-479*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(15,8)=(121*dq[0]*dq[1])/1440;
      upperMatDiffusionTimesMu(15,9)=(121*dq[0]*dq[1])/2880;
      upperMatDiffusionTimesMu(15,10)=(121*dq[0]*dq[1])/1440;
      upperMatDiffusionTimesMu(15,11)=(121*dq[0]*dq[1])/2880;
      upperMatDiffusionTimesMu(15,12)=(-571*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(15,13)=(29*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(15,14)=(-179*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(15,15)=(29*dq[0]*dq[1])/720;
      upperMatDiffusionTimesMu(15,16)=(29*dq[0]*dq[1])/1440;
      upperMatDiffusionTimesMu(15,17)=(-571*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(15,18)=(29*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(15,19)=(-179*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(16,0)=(-479*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(16,1)=(179*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(16,2)=(-1021*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(16,3)=(179*dq[0]*dq[1])/1440;
      upperMatDiffusionTimesMu(16,4)=(179*dq[0]*dq[1])/720;
      upperMatDiffusionTimesMu(16,5)=(-479*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(16,6)=(179*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(16,7)=(-1021*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(16,8)=(121*dq[0]*dq[1])/2880;
      upperMatDiffusionTimesMu(16,9)=(121*dq[0]*dq[1])/1440;
      upperMatDiffusionTimesMu(16,10)=(121*dq[0]*dq[1])/2880;
      upperMatDiffusionTimesMu(16,11)=(121*dq[0]*dq[1])/1440;
      upperMatDiffusionTimesMu(16,12)=(-179*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(16,13)=(29*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(16,14)=(-571*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(16,15)=(29*dq[0]*dq[1])/1440;
      upperMatDiffusionTimesMu(16,16)=(29*dq[0]*dq[1])/720;
      upperMatDiffusionTimesMu(16,17)=(-179*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(16,18)=(29*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(16,19)=(-571*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(17,0)=(101*dq[0]*dq[1])/2304;
      upperMatDiffusionTimesMu(17,1)=(-503*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(17,2)=(1129*dq[0]*dq[1])/23040;
      upperMatDiffusionTimesMu(17,3)=(-1117*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(17,4)=(-503*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(17,5)=(73*dq[0]*dq[1])/1280;
      upperMatDiffusionTimesMu(17,6)=(-1117*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(17,7)=(101*dq[0]*dq[1])/2304;
      upperMatDiffusionTimesMu(17,8)=(-181*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(17,9)=(-151*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(17,10)=-(dq[0]*dq[1])/48;
      upperMatDiffusionTimesMu(17,11)=(-181*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(17,12)=(5*dq[0]*dq[1])/256;
      upperMatDiffusionTimesMu(17,13)=(-73*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(17,14)=(419*dq[0]*dq[1])/23040;
      upperMatDiffusionTimesMu(17,15)=(-49*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(17,16)=(-73*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(17,17)=(23*dq[0]*dq[1])/1280;
      upperMatDiffusionTimesMu(17,18)=(-49*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(17,19)=(5*dq[0]*dq[1])/256;
      upperMatDiffusionTimesMu(18,0)=(-479*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(18,1)=(179*dq[0]*dq[1])/1440;
      upperMatDiffusionTimesMu(18,2)=(-479*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(18,3)=(179*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(18,4)=(179*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(18,5)=(-1021*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(18,6)=(179*dq[0]*dq[1])/720;
      upperMatDiffusionTimesMu(18,7)=(-1021*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(18,8)=(121*dq[0]*dq[1])/2880;
      upperMatDiffusionTimesMu(18,9)=(121*dq[0]*dq[1])/2880;
      upperMatDiffusionTimesMu(18,10)=(121*dq[0]*dq[1])/1440;
      upperMatDiffusionTimesMu(18,11)=(121*dq[0]*dq[1])/1440;
      upperMatDiffusionTimesMu(18,12)=(-179*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(18,13)=(29*dq[0]*dq[1])/1440;
      upperMatDiffusionTimesMu(18,14)=(-179*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(18,15)=(29*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(18,16)=(29*dq[0]*dq[1])/1152;
      upperMatDiffusionTimesMu(18,17)=(-571*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(18,18)=(29*dq[0]*dq[1])/720;
      upperMatDiffusionTimesMu(18,19)=(-571*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(19,0)=(1129*dq[0]*dq[1])/23040;
      upperMatDiffusionTimesMu(19,1)=(-503*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(19,2)=(101*dq[0]*dq[1])/2304;
      upperMatDiffusionTimesMu(19,3)=(-503*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(19,4)=(-1117*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(19,5)=(101*dq[0]*dq[1])/2304;
      upperMatDiffusionTimesMu(19,6)=(-1117*dq[0]*dq[1])/11520;
      upperMatDiffusionTimesMu(19,7)=(73*dq[0]*dq[1])/1280;
      upperMatDiffusionTimesMu(19,8)=(-151*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(19,9)=(-181*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(19,10)=(-181*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(19,11)=-(dq[0]*dq[1])/48;
      upperMatDiffusionTimesMu(19,12)=(419*dq[0]*dq[1])/23040;
      upperMatDiffusionTimesMu(19,13)=(-73*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(19,14)=(5*dq[0]*dq[1])/256;
      upperMatDiffusionTimesMu(19,15)=(-73*dq[0]*dq[1])/5760;
      upperMatDiffusionTimesMu(19,16)=(-49*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(19,17)=(5*dq[0]*dq[1])/256;
      upperMatDiffusionTimesMu(19,18)=(-49*dq[0]*dq[1])/3840;
      upperMatDiffusionTimesMu(19,19)=(23*dq[0]*dq[1])/1280;
    }

    // Multiply relevant matrices by massMatrixInv
    iMatDiffusionTimesMu = massMatrixInv*iMatDiffusionTimesMu;
    lowerMatDiffusionTimesMu = massMatrixInv*lowerMatDiffusionTimesMu;
    upperMatDiffusionTimesMu = massMatrixInv*upperMatDiffusionTimesMu;
  }

  Lucee::UpdaterStatus
  LenardBernsteinDiff3DUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    const Lucee::Field<3, double>& inpFld = this->getInp<Lucee::Field<3, double> >(0);
    const Lucee::Field<3, double>& vThermSqFld = this->getInp<Lucee::Field<3, double> >(1);
    Lucee::Field<3, double>& diffOut = this->getOut<Lucee::Field<3, double> >(0);

    double dt = t-this->getCurrTime();
    
    Lucee::ConstFieldPtr<double> inpFldPtr  = inpFld.createConstPtr();
    Lucee::ConstFieldPtr<double> vThermSqPtr  = vThermSqFld.createConstPtr();
    Lucee::FieldPtr<double> diffOutPtr = diffOut.createPtr();

    // check time-step
    double cflm = 1.1*cfl;
    double dxMax2 = grid.getDx(1)*grid.getDx(1);
    double cfla = 0.0;
    
    diffOut = 0.0;

    // local region to index
    Lucee::Region<3, int> localRgn = grid.getLocalRegion();
    Lucee::Region<3, int> globalRgn = grid.getGlobalRegion();

    Lucee::RowMajorSequencer<3> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[3];
    double xc[3];
    seq.fillWithIndex(idx);
    nodalBasis->setIndex(idx);
    unsigned nlocal = nodalBasis->getNumNodes(); 

    // Get alpha
    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;
    std::vector<double> resultVector(1);
    evaluateFunction(*L, t, resultVector);
    alpha = resultVector[0];

    for (int ix = localRgn.getLower(0); ix < localRgn.getUpper(0); ix++)
    {
      for (int iv = localRgn.getLower(1); iv < localRgn.getUpper(1); iv++)
      {
        for (int iMu = localRgn.getLower(2); iMu < localRgn.getUpper(2); iMu++)
        {
          diffOut.setPtr(diffOutPtr, ix, iv, iMu);
          inpFld.setPtr(inpFldPtr, ix, iv, iMu);
          vThermSqFld.setPtr(vThermSqPtr, ix, iv, iMu);
          
          // Keep track of maximum cfla (v-parallel diffusion)
          cfla = std::max(cfla, alpha*vThermSqPtr[0]*dt/(grid.getDx(1)*grid.getDx(1)) );
          // Keep track of maximum cfla (mu diffusion)
          double muUpper = xc[2] + grid.getDx(2);
          cfla = std::max(cfla, alpha*speciesMass*vThermSqPtr[0]/B0*dt/muUpper*
              std::max(1.0, 4*muUpper*muUpper/(grid.getDx(2)*grid.getDx(2))) );
          
          if (cfla > cflm)
            return Lucee::UpdaterStatus(false, dt*cfl/cfla);
          
          grid.setIndex(ix, iv, iMu);
          grid.getCentroid(xc);


          Eigen::VectorXd fAtNodes(nlocal);

          for (int i = 0; i < nlocal; i++)
            fAtNodes(i) = inpFldPtr[i];

          //Eigen::VectorXd updateF = 2*speciesMass/B0*(xc[2]*iMat[2] + iMatDiffusionTimesMu)*fAtNodes;

          Eigen::VectorXd updateF = ( iMat[1] + 
              2*speciesMass/B0*(xc[2]*iMat[2] + iMatDiffusionTimesMu) )*fAtNodes;
          
          Eigen::VectorXd fLowerAtNodes(nlocal);
          Eigen::VectorXd fUpperAtNodes(nlocal);

          // add in contribution from cells attached to lower/upper faces in vPara
          inpFld.setPtr(inpFldPtr, ix, iv-1, iMu); // cell attached to lower face
          for (int i = 0; i < nlocal; i++)
            fLowerAtNodes(i) = inpFldPtr[i];
          updateF = updateF + lowerMat[1]*fLowerAtNodes;

          inpFld.setPtr(inpFldPtr, ix, iv+1, iMu); // cell attached to upper face
          for (int i = 0; i < nlocal; i++)
            fUpperAtNodes(i) = inpFldPtr[i];
          updateF = updateF + upperMat[1]*fUpperAtNodes;
          
          // add in contribution from cells attached to lower/upper faces in mu
          inpFld.setPtr(inpFldPtr, ix, iv, iMu-1); // cell attached to lower face
          for (int i = 0; i < nlocal; i++)
            fLowerAtNodes(i) = inpFldPtr[i];
          updateF = updateF + 2*speciesMass/B0*(xc[2]*lowerMat[2] + lowerMatDiffusionTimesMu)*fLowerAtNodes;

          inpFld.setPtr(inpFldPtr, ix, iv, iMu+1); // cell attached to upper face
          for (int i = 0; i < nlocal; i++)
            fUpperAtNodes(i) = inpFldPtr[i];
          updateF = updateF + 2*speciesMass/B0*(xc[2]*upperMat[2] + upperMatDiffusionTimesMu)*fUpperAtNodes;

          // Accumulate updateF to output
          for (int i = 0; i < nlocal; i++)
            diffOutPtr[i] = diffOutPtr[i] + updateF(i);
        }
      }
    }

    seq.reset();//= Lucee::RowMajorSequencer<3>(localRgn);
    // Final sweep, update solution with forward Euler step

    while (seq.step())
    {
      seq.fillWithIndex(idx);
      diffOut.setPtr(diffOutPtr, idx);
      vThermSqFld.setPtr(vThermSqPtr, idx);
      
      if (onlyIncrement == false)
      {
        inpFld.setPtr(inpFldPtr, idx);
        for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
          diffOutPtr[nodeIndex] = inpFldPtr[nodeIndex] +
            dt*alpha*vThermSqPtr[nodeIndex]*diffOutPtr[nodeIndex];
      }
      else
      {
        for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
          diffOutPtr[nodeIndex] = alpha*vThermSqPtr[nodeIndex]*diffOutPtr[nodeIndex];
      }
    }

    return Lucee::UpdaterStatus(true, dt*cfl/cfla);
  }

  void
  LenardBernsteinDiff3DUpdater::declareTypes()
  {
    // Input: distribution function
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: vThermSq field
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Output: distribution function
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
  }

  void
  LenardBernsteinDiff3DUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
    {
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
      {
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
      }
    }
  }

  void
  LenardBernsteinDiff3DUpdater::evaluateFunction(Lucee::LuaState& L, double tm,
    std::vector<double>& res)
  {
    // push function object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, fnRef);
    // push variables on stack
    lua_pushnumber(L, tm);
    // call function
    if (lua_pcall(L, 1, res.size(), 0) != 0)
    {
      Lucee::Except lce("HeatFluxAtEdgeUpdater::evaluateFunction: ");
      lce << "Problem evaluating function supplied as 'tPerpProfile' "
          << std::endl;
      std::string err(lua_tostring(L, -1));
      lua_pop(L, 1);
      lce << "[" << err << "]";
      throw lce;
    }
    // fetch results
    for (int i=-res.size(); i<0; ++i)
    {
      if (!lua_isnumber(L, i))
        throw Lucee::Except("HeatFluxAtEdgeUpdater::evaluateFunction: Return value not a number");
      res[res.size()+i] = lua_tonumber(L, i);
    }
    lua_pop(L, 1);
  }
}
