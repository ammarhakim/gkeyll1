/**
 * @file	LcLenardBernsteinDiffPara3D2VUpdater.cpp
 *
 * @brief	Updater to evaluate the diffusion term in the L-B collision operator for 3D2V problems.
 * This updater is different from LcLenardBernsteinDiff5DUpdater in its use of a 2D recovery method
 * to evaluate diffusion terms
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLenardBernsteinDiffPara3D2VUpdater.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  using namespace Eigen;
  const char *LenardBernsteinDiffPara3D2VUpdater::id = "LenardBernsteinDiffPara3D2VUpdater";

  LenardBernsteinDiffPara3D2VUpdater::LenardBernsteinDiffPara3D2VUpdater()
  {
  }

  LenardBernsteinDiffPara3D2VUpdater::~LenardBernsteinDiffPara3D2VUpdater()
  {
  }

  void
  LenardBernsteinDiffPara3D2VUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except("LenardBernsteinDiffPara3D2VUpdater::readInput: Must specify element to use using 'basis5d'");
 
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("LenardBernsteinDiffPara3D2VUpdater::readInput: Must specify element to use using 'basis3d'");
 

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis2d"))
      nodalBasis2d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis2d");
    else
      throw Lucee::Except("LenardBernsteinDiffPara3D2VUpdater::readInput: Must specify element to use using 'basis2d'");
 
    // CFL number to control time-step
    cfl = tbl.getNumber("cfl"); // CFL number

    onlyIncrement = false;
    if (tbl.hasBool("onlyIncrement"))
      onlyIncrement = tbl.getBool("onlyIncrement");

    if (tbl.hasNumber("speciesMass"))
      speciesMass = tbl.getNumber("speciesMass");
    else
      throw Lucee::Except("LenardBernsteinDiffPara3D2VUpdater::readInput: Must specify speciesMass");

    if (tbl.hasFunction("alpha"))
      fnRef = tbl.getFunctionRef("alpha");
    else
      throw Lucee::Except("LenardBernsteinDiffPara3D2VUpdater::readInput: Must supply a collision frequency function as alpha.");
  }

  void
  LenardBernsteinDiffPara3D2VUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();
    
// get hold of grid
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();
    // local region to update
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();

    Lucee::RowMajorSequencer<5> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[5];
    seq.fillWithIndex(idx);
    nodalBasis5d->setIndex(idx);
    grid.setIndex(idx);
    nodalBasis2d->setIndex(idx[0],idx[1]);

    unsigned nlocal5d = nodalBasis5d->getNumNodes();
    unsigned nlocal2d = nodalBasis2d->getNumNodes();

    // Get a copy of the nodal coordinates
    Lucee::Matrix<double> nodeCoordsLucee(nlocal5d, 5);
    nodalBasis5d->getNodalCoordinates(nodeCoordsLucee);
    Eigen::MatrixXd nodeCoords(nlocal5d, 5);
    copyLuceeToEigen(nodeCoordsLucee, nodeCoords);

    double dxMin = grid.getDx(0);
    for (int d = 1; d < 3; d++)
      dxMin = std::min(dxMin, grid.getDx(d));

    // Find all nodes that share the same location as node zero
    // Will eventually need to do this at all nodes on lower surface
    // but can get away with doing this at one point for linear element
    nodalStencil = std::vector<int>(nlocal5d);
    int stencilIndex = 0;
    for (int nodeIndex = 0; nodeIndex < nlocal5d; nodeIndex++)
    {
      if (sameConfigCoords(0, nodeIndex, dxMin, nodeCoords) == true)
      {
        nodalStencil[stencilIndex] = nodeIndex;
        stencilIndex++;
      }
    }
    nodalStencil.resize(stencilIndex);

    lowerMat.resize(2, Eigen::MatrixXd(nlocal2d, nlocal2d));
    upperMat.resize(2, Eigen::MatrixXd(nlocal2d, nlocal2d));
    upperCenter.resize(2, Eigen::MatrixXd(nlocal2d, nlocal2d));
    selfCenter.resize(2, Eigen::MatrixXd(nlocal2d, nlocal2d));
    lowerCenter.resize(2, Eigen::MatrixXd(nlocal2d, nlocal2d));

    upperCenterTimesMu = Eigen::MatrixXd(nlocal2d, nlocal2d);
    selfCenterTimesMu = Eigen::MatrixXd(nlocal2d, nlocal2d);
    lowerCenterTimesMu = Eigen::MatrixXd(nlocal2d, nlocal2d);
    lowerMatTimesMu = Eigen::MatrixXd(nlocal2d, nlocal2d);
    upperMatTimesMu = Eigen::MatrixXd(nlocal2d, nlocal2d);

    double dq[2];
    dq[0] = grid.getDx(3);
    dq[1] = grid.getDx(4);

    // Code generated from mathematica
    upperCenter[0](0,0)=dq[1]/(36*dq[0]);
    upperCenter[0](0,1)=(5*dq[1])/(36*dq[0]);
    upperCenter[0](0,2)=dq[1]/(72*dq[0]);
    upperCenter[0](0,3)=(5*dq[1])/(72*dq[0]);
    upperCenter[0](1,0)=(-7*dq[1])/(36*dq[0]);
    upperCenter[0](1,1)=(-13*dq[1])/(18*dq[0]);
    upperCenter[0](1,2)=(-7*dq[1])/(72*dq[0]);
    upperCenter[0](1,3)=(-13*dq[1])/(36*dq[0]);
    upperCenter[0](2,0)=dq[1]/(72*dq[0]);
    upperCenter[0](2,1)=(5*dq[1])/(72*dq[0]);
    upperCenter[0](2,2)=dq[1]/(36*dq[0]);
    upperCenter[0](2,3)=(5*dq[1])/(36*dq[0]);
    upperCenter[0](3,0)=(-7*dq[1])/(72*dq[0]);
    upperCenter[0](3,1)=(-13*dq[1])/(36*dq[0]);
    upperCenter[0](3,2)=(-7*dq[1])/(36*dq[0]);
    upperCenter[0](3,3)=(-13*dq[1])/(18*dq[0]);
    upperMat[0](0,0)=(5*dq[1])/(36*dq[0]);
    upperMat[0](0,1)=dq[1]/(36*dq[0]);
    upperMat[0](0,2)=(5*dq[1])/(72*dq[0]);
    upperMat[0](0,3)=dq[1]/(72*dq[0]);
    upperMat[0](1,0)=(4*dq[1])/(9*dq[0]);
    upperMat[0](1,1)=(5*dq[1])/(36*dq[0]);
    upperMat[0](1,2)=(2*dq[1])/(9*dq[0]);
    upperMat[0](1,3)=(5*dq[1])/(72*dq[0]);
    upperMat[0](2,0)=(5*dq[1])/(72*dq[0]);
    upperMat[0](2,1)=dq[1]/(72*dq[0]);
    upperMat[0](2,2)=(5*dq[1])/(36*dq[0]);
    upperMat[0](2,3)=dq[1]/(36*dq[0]);
    upperMat[0](3,0)=(2*dq[1])/(9*dq[0]);
    upperMat[0](3,1)=(5*dq[1])/(72*dq[0]);
    upperMat[0](3,2)=(4*dq[1])/(9*dq[0]);
    upperMat[0](3,3)=(5*dq[1])/(36*dq[0]);
    selfCenter[0](0,0)=0;
    selfCenter[0](0,1)=0;
    selfCenter[0](0,2)=0;
    selfCenter[0](0,3)=0;
    selfCenter[0](1,0)=0;
    selfCenter[0](1,1)=0;
    selfCenter[0](1,2)=0;
    selfCenter[0](1,3)=0;
    selfCenter[0](2,0)=0;
    selfCenter[0](2,1)=0;
    selfCenter[0](2,2)=0;
    selfCenter[0](2,3)=0;
    selfCenter[0](3,0)=0;
    selfCenter[0](3,1)=0;
    selfCenter[0](3,2)=0;
    selfCenter[0](3,3)=0;
    lowerCenter[0](0,0)=(-13*dq[1])/(18*dq[0]);
    lowerCenter[0](0,1)=(-7*dq[1])/(36*dq[0]);
    lowerCenter[0](0,2)=(-13*dq[1])/(36*dq[0]);
    lowerCenter[0](0,3)=(-7*dq[1])/(72*dq[0]);
    lowerCenter[0](1,0)=(5*dq[1])/(36*dq[0]);
    lowerCenter[0](1,1)=dq[1]/(36*dq[0]);
    lowerCenter[0](1,2)=(5*dq[1])/(72*dq[0]);
    lowerCenter[0](1,3)=dq[1]/(72*dq[0]);
    lowerCenter[0](2,0)=(-13*dq[1])/(36*dq[0]);
    lowerCenter[0](2,1)=(-7*dq[1])/(72*dq[0]);
    lowerCenter[0](2,2)=(-13*dq[1])/(18*dq[0]);
    lowerCenter[0](2,3)=(-7*dq[1])/(36*dq[0]);
    lowerCenter[0](3,0)=(5*dq[1])/(72*dq[0]);
    lowerCenter[0](3,1)=dq[1]/(72*dq[0]);
    lowerCenter[0](3,2)=(5*dq[1])/(36*dq[0]);
    lowerCenter[0](3,3)=dq[1]/(36*dq[0]);
    lowerMat[0](0,0)=(5*dq[1])/(36*dq[0]);
    lowerMat[0](0,1)=(4*dq[1])/(9*dq[0]);
    lowerMat[0](0,2)=(5*dq[1])/(72*dq[0]);
    lowerMat[0](0,3)=(2*dq[1])/(9*dq[0]);
    lowerMat[0](1,0)=dq[1]/(36*dq[0]);
    lowerMat[0](1,1)=(5*dq[1])/(36*dq[0]);
    lowerMat[0](1,2)=dq[1]/(72*dq[0]);
    lowerMat[0](1,3)=(5*dq[1])/(72*dq[0]);
    lowerMat[0](2,0)=(5*dq[1])/(72*dq[0]);
    lowerMat[0](2,1)=(2*dq[1])/(9*dq[0]);
    lowerMat[0](2,2)=(5*dq[1])/(36*dq[0]);
    lowerMat[0](2,3)=(4*dq[1])/(9*dq[0]);
    lowerMat[0](3,0)=dq[1]/(72*dq[0]);
    lowerMat[0](3,1)=(5*dq[1])/(72*dq[0]);
    lowerMat[0](3,2)=dq[1]/(36*dq[0]);
    lowerMat[0](3,3)=(5*dq[1])/(36*dq[0]);
    upperCenter[1](0,0)=dq[0]/(36*dq[1]);
    upperCenter[1](0,1)=dq[0]/(72*dq[1]);
    upperCenter[1](0,2)=(5*dq[0])/(36*dq[1]);
    upperCenter[1](0,3)=(5*dq[0])/(72*dq[1]);
    upperCenter[1](1,0)=dq[0]/(72*dq[1]);
    upperCenter[1](1,1)=dq[0]/(36*dq[1]);
    upperCenter[1](1,2)=(5*dq[0])/(72*dq[1]);
    upperCenter[1](1,3)=(5*dq[0])/(36*dq[1]);
    upperCenter[1](2,0)=(-7*dq[0])/(36*dq[1]);
    upperCenter[1](2,1)=(-7*dq[0])/(72*dq[1]);
    upperCenter[1](2,2)=(-13*dq[0])/(18*dq[1]);
    upperCenter[1](2,3)=(-13*dq[0])/(36*dq[1]);
    upperCenter[1](3,0)=(-7*dq[0])/(72*dq[1]);
    upperCenter[1](3,1)=(-7*dq[0])/(36*dq[1]);
    upperCenter[1](3,2)=(-13*dq[0])/(36*dq[1]);
    upperCenter[1](3,3)=(-13*dq[0])/(18*dq[1]);
    upperMat[1](0,0)=(5*dq[0])/(36*dq[1]);
    upperMat[1](0,1)=(5*dq[0])/(72*dq[1]);
    upperMat[1](0,2)=dq[0]/(36*dq[1]);
    upperMat[1](0,3)=dq[0]/(72*dq[1]);
    upperMat[1](1,0)=(5*dq[0])/(72*dq[1]);
    upperMat[1](1,1)=(5*dq[0])/(36*dq[1]);
    upperMat[1](1,2)=dq[0]/(72*dq[1]);
    upperMat[1](1,3)=dq[0]/(36*dq[1]);
    upperMat[1](2,0)=(4*dq[0])/(9*dq[1]);
    upperMat[1](2,1)=(2*dq[0])/(9*dq[1]);
    upperMat[1](2,2)=(5*dq[0])/(36*dq[1]);
    upperMat[1](2,3)=(5*dq[0])/(72*dq[1]);
    upperMat[1](3,0)=(2*dq[0])/(9*dq[1]);
    upperMat[1](3,1)=(4*dq[0])/(9*dq[1]);
    upperMat[1](3,2)=(5*dq[0])/(72*dq[1]);
    upperMat[1](3,3)=(5*dq[0])/(36*dq[1]);
    lowerCenter[1](0,0)=(-13*dq[0])/(18*dq[1]);
    lowerCenter[1](0,1)=(-13*dq[0])/(36*dq[1]);
    lowerCenter[1](0,2)=(-7*dq[0])/(36*dq[1]);
    lowerCenter[1](0,3)=(-7*dq[0])/(72*dq[1]);
    lowerCenter[1](1,0)=(-13*dq[0])/(36*dq[1]);
    lowerCenter[1](1,1)=(-13*dq[0])/(18*dq[1]);
    lowerCenter[1](1,2)=(-7*dq[0])/(72*dq[1]);
    lowerCenter[1](1,3)=(-7*dq[0])/(36*dq[1]);
    lowerCenter[1](2,0)=(5*dq[0])/(36*dq[1]);
    lowerCenter[1](2,1)=(5*dq[0])/(72*dq[1]);
    lowerCenter[1](2,2)=dq[0]/(36*dq[1]);
    lowerCenter[1](2,3)=dq[0]/(72*dq[1]);
    lowerCenter[1](3,0)=(5*dq[0])/(72*dq[1]);
    lowerCenter[1](3,1)=(5*dq[0])/(36*dq[1]);
    lowerCenter[1](3,2)=dq[0]/(72*dq[1]);
    lowerCenter[1](3,3)=dq[0]/(36*dq[1]);
    lowerMat[1](0,0)=(5*dq[0])/(36*dq[1]);
    lowerMat[1](0,1)=(5*dq[0])/(72*dq[1]);
    lowerMat[1](0,2)=(4*dq[0])/(9*dq[1]);
    lowerMat[1](0,3)=(2*dq[0])/(9*dq[1]);
    lowerMat[1](1,0)=(5*dq[0])/(72*dq[1]);
    lowerMat[1](1,1)=(5*dq[0])/(36*dq[1]);
    lowerMat[1](1,2)=(2*dq[0])/(9*dq[1]);
    lowerMat[1](1,3)=(4*dq[0])/(9*dq[1]);
    lowerMat[1](2,0)=dq[0]/(36*dq[1]);
    lowerMat[1](2,1)=dq[0]/(72*dq[1]);
    lowerMat[1](2,2)=(5*dq[0])/(36*dq[1]);
    lowerMat[1](2,3)=(5*dq[0])/(72*dq[1]);
    lowerMat[1](3,0)=dq[0]/(72*dq[1]);
    lowerMat[1](3,1)=dq[0]/(36*dq[1]);
    lowerMat[1](3,2)=(5*dq[0])/(72*dq[1]);
    lowerMat[1](3,3)=(5*dq[0])/(36*dq[1]);
    selfCenter[1](0,0)=0;
    selfCenter[1](0,1)=0;
    selfCenter[1](0,2)=0;
    selfCenter[1](0,3)=0;
    selfCenter[1](1,0)=0;
    selfCenter[1](1,1)=0;
    selfCenter[1](1,2)=0;
    selfCenter[1](1,3)=0;
    selfCenter[1](2,0)=0;
    selfCenter[1](2,1)=0;
    selfCenter[1](2,2)=0;
    selfCenter[1](2,3)=0;
    selfCenter[1](3,0)=0;
    selfCenter[1](3,1)=0;
    selfCenter[1](3,2)=0;
    selfCenter[1](3,3)=0;
    upperCenterTimesMu(0,0)=dq[0]/72;
    upperCenterTimesMu(0,1)=dq[0]/144;
    upperCenterTimesMu(0,2)=(5*dq[0])/72;
    upperCenterTimesMu(0,3)=(5*dq[0])/144;
    upperCenterTimesMu(1,0)=dq[0]/144;
    upperCenterTimesMu(1,1)=dq[0]/72;
    upperCenterTimesMu(1,2)=(5*dq[0])/144;
    upperCenterTimesMu(1,3)=(5*dq[0])/72;
    upperCenterTimesMu(2,0)=(-7*dq[0])/72;
    upperCenterTimesMu(2,1)=(-7*dq[0])/144;
    upperCenterTimesMu(2,2)=(-13*dq[0])/36;
    upperCenterTimesMu(2,3)=(-13*dq[0])/72;
    upperCenterTimesMu(3,0)=(-7*dq[0])/144;
    upperCenterTimesMu(3,1)=(-7*dq[0])/72;
    upperCenterTimesMu(3,2)=(-13*dq[0])/72;
    upperCenterTimesMu(3,3)=(-13*dq[0])/36;
    upperMatTimesMu(0,0)=(5*dq[0])/72;
    upperMatTimesMu(0,1)=(5*dq[0])/144;
    upperMatTimesMu(0,2)=dq[0]/72;
    upperMatTimesMu(0,3)=dq[0]/144;
    upperMatTimesMu(1,0)=(5*dq[0])/144;
    upperMatTimesMu(1,1)=(5*dq[0])/72;
    upperMatTimesMu(1,2)=dq[0]/144;
    upperMatTimesMu(1,3)=dq[0]/72;
    upperMatTimesMu(2,0)=(2*dq[0])/9;
    upperMatTimesMu(2,1)=dq[0]/9;
    upperMatTimesMu(2,2)=(5*dq[0])/72;
    upperMatTimesMu(2,3)=(5*dq[0])/144;
    upperMatTimesMu(3,0)=dq[0]/9;
    upperMatTimesMu(3,1)=(2*dq[0])/9;
    upperMatTimesMu(3,2)=(5*dq[0])/144;
    upperMatTimesMu(3,3)=(5*dq[0])/72;
    lowerCenterTimesMu(0,0)=(13*dq[0])/36;
    lowerCenterTimesMu(0,1)=(13*dq[0])/72;
    lowerCenterTimesMu(0,2)=(7*dq[0])/72;
    lowerCenterTimesMu(0,3)=(7*dq[0])/144;
    lowerCenterTimesMu(1,0)=(13*dq[0])/72;
    lowerCenterTimesMu(1,1)=(13*dq[0])/36;
    lowerCenterTimesMu(1,2)=(7*dq[0])/144;
    lowerCenterTimesMu(1,3)=(7*dq[0])/72;
    lowerCenterTimesMu(2,0)=(-5*dq[0])/72;
    lowerCenterTimesMu(2,1)=(-5*dq[0])/144;
    lowerCenterTimesMu(2,2)=-dq[0]/72;
    lowerCenterTimesMu(2,3)=-dq[0]/144;
    lowerCenterTimesMu(3,0)=(-5*dq[0])/144;
    lowerCenterTimesMu(3,1)=(-5*dq[0])/72;
    lowerCenterTimesMu(3,2)=-dq[0]/144;
    lowerCenterTimesMu(3,3)=-dq[0]/72;
    lowerMatTimesMu(0,0)=(-5*dq[0])/72;
    lowerMatTimesMu(0,1)=(-5*dq[0])/144;
    lowerMatTimesMu(0,2)=(-2*dq[0])/9;
    lowerMatTimesMu(0,3)=-dq[0]/9;
    lowerMatTimesMu(1,0)=(-5*dq[0])/144;
    lowerMatTimesMu(1,1)=(-5*dq[0])/72;
    lowerMatTimesMu(1,2)=-dq[0]/9;
    lowerMatTimesMu(1,3)=(-2*dq[0])/9;
    lowerMatTimesMu(2,0)=-dq[0]/72;
    lowerMatTimesMu(2,1)=-dq[0]/144;
    lowerMatTimesMu(2,2)=(-5*dq[0])/72;
    lowerMatTimesMu(2,3)=(-5*dq[0])/144;
    lowerMatTimesMu(3,0)=-dq[0]/144;
    lowerMatTimesMu(3,1)=-dq[0]/72;
    lowerMatTimesMu(3,2)=(-5*dq[0])/144;
    lowerMatTimesMu(3,3)=(-5*dq[0])/72;
    selfCenterTimesMu(0,0)=-dq[0]/6;
    selfCenterTimesMu(0,1)=-dq[0]/12;
    selfCenterTimesMu(0,2)=-dq[0]/6;
    selfCenterTimesMu(0,3)=-dq[0]/12;
    selfCenterTimesMu(1,0)=-dq[0]/12;
    selfCenterTimesMu(1,1)=-dq[0]/6;
    selfCenterTimesMu(1,2)=-dq[0]/12;
    selfCenterTimesMu(1,3)=-dq[0]/6;
    selfCenterTimesMu(2,0)=dq[0]/6;
    selfCenterTimesMu(2,1)=dq[0]/12;
    selfCenterTimesMu(2,2)=dq[0]/6;
    selfCenterTimesMu(2,3)=dq[0]/12;
    selfCenterTimesMu(3,0)=dq[0]/12;
    selfCenterTimesMu(3,1)=dq[0]/6;
    selfCenterTimesMu(3,2)=dq[0]/12;
    selfCenterTimesMu(3,3)=dq[0]/6;

    // pre-multiply each of the matrices by inverse matrix
    Lucee::Matrix<double> massMatrixLucee(nlocal2d, nlocal2d);
    nodalBasis2d->getMassMatrix(massMatrixLucee);
    Eigen::MatrixXd massMatrix(nlocal2d, nlocal2d);
    copyLuceeToEigen(massMatrixLucee, massMatrix);

    // Scale massMatrix so that it can be used for the (v,mu) grid
    massMatrix *= grid.getDx(3)*grid.getDx(4)/(grid.getDx(0)*grid.getDx(1));

    Eigen::MatrixXd massMatrixInv = massMatrix.inverse();

    lowerMat[0] = massMatrixInv*lowerMat[0];
    upperMat[0] = massMatrixInv*upperMat[0];
    upperCenter[0] = massMatrixInv*upperCenter[0];
    selfCenter[0] = massMatrixInv*selfCenter[0];
    lowerCenter[0] = massMatrixInv*lowerCenter[0];
    
    lowerMat[1] = massMatrixInv*lowerMat[1];
    upperMat[1] = massMatrixInv*upperMat[1];
    upperCenter[1] = massMatrixInv*upperCenter[1];
    selfCenter[1] = massMatrixInv*selfCenter[1];
    lowerCenter[1] = massMatrixInv*lowerCenter[1];

    lowerMatTimesMu = massMatrixInv*lowerMatTimesMu;
    upperMatTimesMu = massMatrixInv*upperMatTimesMu;
    upperCenterTimesMu = massMatrixInv*upperCenterTimesMu;
    selfCenterTimesMu = massMatrixInv*selfCenterTimesMu;
    lowerCenterTimesMu = massMatrixInv*lowerCenterTimesMu;
  }

  Lucee::UpdaterStatus
  LenardBernsteinDiffPara3D2VUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // Distribution function
    const Lucee::Field<5, double>& fIn = this->getInp<Lucee::Field<5, double> >(0);
    // Temperature in joules
    const Lucee::Field<3, double>& temperatureIn = this->getInp<Lucee::Field<3, double> >(1);
    // Magnetic field
    const Lucee::Field<3, double>& bFieldIn = this->getInp<Lucee::Field<3, double> >(2);
    // Dimensionally correct number density from weighted moment calculation
    const Lucee::Field<3, double>& numDensityIn = this->getInp<Lucee::Field<3, double> >(3);
    // Output distribution function
    Lucee::Field<5, double>& fOut = this->getOut<Lucee::Field<5, double> >(0);

    double dt = t-this->getCurrTime();

    // Input fields
    Lucee::ConstFieldPtr<double> fInPtr = fIn.createConstPtr();
    Lucee::ConstFieldPtr<double> temperatureInPtr = temperatureIn.createConstPtr();
    Lucee::ConstFieldPtr<double> bFieldInPtr = bFieldIn.createConstPtr();
    Lucee::ConstFieldPtr<double> numDensityInPtr = numDensityIn.createConstPtr();
    // Writeable fields
    Lucee::FieldPtr<double> fOutPtr = fOut.createPtr();

    // check time-step
    double cflm = 1.1*cfl;
    double cfla = 0.0;
    
    fOut = 0.0;

    // local region to index
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();
    Lucee::Region<5, int> globalRgn = grid.getGlobalRegion();

    Lucee::RowMajorSequencer<5> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[5];
    double cellCentroid[5];
    seq.fillWithIndex(idx);
    nodalBasis5d->setIndex(idx);
    int nlocal5d = nodalBasis5d->getNumNodes(); 
    nodalBasis3d->setIndex(idx[0], idx[1], idx[2]);
    int nlocal3d = nodalBasis3d->getNumNodes(); 
    nodalBasis2d->setIndex(idx[0], idx[1]);
    int nlocal2d = nodalBasis2d->getNumNodes(); 

    // Get alpha. Need to scale by n/T^(3/2) in this class
    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;
    std::vector<double> resultVector(1);
    evaluateFunction(*L, t, resultVector);
    double alpha = resultVector[0];

    // Should be size 4 for linear elements
    Eigen::VectorXd fReduced(nodalStencil.size());
    Eigen::VectorXd fLowerReduced(nodalStencil.size());
    Eigen::VectorXd fUpperReduced(nodalStencil.size());

    for (int ix = localRgn.getLower(0); ix < localRgn.getUpper(0); ix++)
    {
      idx[0] = ix;
      for (int iy = localRgn.getLower(1); iy < localRgn.getUpper(1); iy++)
      {
        idx[1] = iy;
        for (int iz = localRgn.getLower(2); iz < localRgn.getUpper(2); iz++)
        {
          idx[2] = iz;
          bFieldIn.setPtr(bFieldInPtr, idx[0], idx[1], idx[2]);
          temperatureIn.setPtr(temperatureInPtr, idx[0], idx[1], idx[2]);
          numDensityIn.setPtr(numDensityInPtr, idx[0], idx[1], idx[2]);

          // At this location, loop over each configuration space grid node
          for (int configNode = 0; configNode < nlocal3d; configNode++)
          {
            // Keep track of max CFL number
            cfla = std::max( cfla, std::abs(4.0*alpha*numDensityInPtr[configNode]/(temperatureInPtr[configNode]*sqrt(temperatureInPtr[configNode]))*
              temperatureInPtr[configNode]/speciesMass*dt/(grid.getDx(3)*grid.getDx(3))) );

            // Loop over entire (vPar,mu) space
            for (int iv = localRgn.getLower(3); iv < localRgn.getUpper(3); iv++)
            {
              idx[3] = iv;
              for (int iMu = localRgn.getLower(4); iMu < localRgn.getUpper(4); iMu++)
              {
                idx[4] = iMu;

                // Set pointers to current location on grid
                fIn.setPtr(fInPtr, idx);
                fOut.setPtr(fOutPtr, idx);
                // Get the coordinates of cell center
                grid.setIndex(idx);
                grid.getCentroid(cellCentroid);

                // Fill out fReduced at this location
                for (int i = 0; i < fReduced.size(); i++)
                  fReduced(i) = fInPtr[configNode + nodalStencil[i]];

                Eigen::VectorXd updateF = selfCenter[0]*fReduced;

                // add in contribution from cells attached to lower/upper faces in vPara
                if (iv > globalRgn.getLower(3))
                {
                  updateF = updateF + lowerCenter[0]*fReduced;
                  idx[3] = idx[3] - 1;
                  fIn.setPtr(fInPtr, idx); // cell attached to lower face
                  for (int i = 0; i < fLowerReduced.size(); i++)
                    fLowerReduced(i) = fInPtr[configNode + nodalStencil[i]];
                  updateF = updateF + lowerMat[0]*fLowerReduced;
                  idx[3] = idx[3] + 1;
                }

                if (iv < globalRgn.getUpper(3)-1)
                {
                  updateF = updateF + upperCenter[0]*fReduced;
                  idx[3] = idx[3] + 1;
                  fIn.setPtr(fInPtr, idx); // cell attached to upper face
                  for (int i = 0; i < fUpperReduced.size(); i++)
                    fUpperReduced(i) = fInPtr[configNode + nodalStencil[i]];
                  updateF = updateF + upperMat[0]*fUpperReduced;
                  idx[3] = idx[3] - 1;
                }

                // Accumulate updateF to output
                for (int i = 0; i < fReduced.size(); i++)
                  fOutPtr[configNode + nodalStencil[i]] = fOutPtr[configNode + nodalStencil[i]] + updateF(i);
              }
            }
          }
        }
      }
    }
    
    if (cfla > cflm)
      return Lucee::UpdaterStatus(false, dt*cfl/cfla);

    seq.reset();
    // Final sweep, update solution with forward Euler step
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      fOut.setPtr(fOutPtr, idx);
      // Compute alpha scale factor n/Te^(3/2)
      temperatureIn.setPtr(temperatureInPtr, idx[0], idx[1], idx[2]);
      numDensityIn.setPtr(numDensityInPtr, idx[0], idx[1], idx[2]);

      if (onlyIncrement == false)
      {
        fIn.setPtr(fInPtr, idx);
        for (int configNode = 0; configNode < nlocal3d; configNode++)
        {
          for (int stencilIndex = 0; stencilIndex < nodalStencil.size(); stencilIndex++)
          {
            fOutPtr[configNode + nodalStencil[stencilIndex]] = fInPtr[configNode + nodalStencil[stencilIndex]]
              + dt*alpha*numDensityInPtr[configNode]/(temperatureInPtr[configNode]*sqrt(temperatureInPtr[configNode]))*
                  temperatureInPtr[configNode]/speciesMass*fOutPtr[configNode + nodalStencil[stencilIndex]];
          }
        }
      }
      else
      {
        for (int configNode = 0; configNode < nlocal3d; configNode++)
        {
          for (int stencilIndex = 0; stencilIndex < nodalStencil.size(); stencilIndex++)
          {
            fOutPtr[configNode + nodalStencil[stencilIndex]] = alpha*numDensityInPtr[configNode]/(temperatureInPtr[configNode]
                *sqrt(temperatureInPtr[configNode]))
                *temperatureInPtr[configNode]/speciesMass*fOutPtr[configNode + nodalStencil[stencilIndex]];
          }
        }
      }
    }

    return Lucee::UpdaterStatus(true, dt*cfl/cfla);
  }

  void
  LenardBernsteinDiffPara3D2VUpdater::declareTypes()
  {
    // Input: Distribution function
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Input: Temperature in joules
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: Magnetic field
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: dimensionally correct number density
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // returns one output (fNew)
    this->appendOutVarType(typeid(Lucee::Field<5, double>));
  }

  void
  LenardBernsteinDiffPara3D2VUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
  LenardBernsteinDiffPara3D2VUpdater::evaluateFunction(Lucee::LuaState& L, double tm,
    std::vector<double>& res)
  {
    // push function object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, fnRef);
    // push variables on stack
    lua_pushnumber(L, tm);
    // call function
    if (lua_pcall(L, 1, res.size(), 0) != 0)
    {
      Lucee::Except lce("LenardBernsteinDiffPara3D2VUpdater::evaluateFunction: ");
      lce << "Problem evaluating function supplied as 'alpha' "
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
        throw Lucee::Except("LenardBernsteinDiffPara3D2VUpdater::evaluateFunction: Return value not a number");
      res[res.size()+i] = lua_tonumber(L, i);
    }
    lua_pop(L, 1);
  }

  bool
  LenardBernsteinDiffPara3D2VUpdater::sameConfigCoords(int srcIndex, int tarIndex, double dxMin,
    const Eigen::MatrixXd& nodeList)
  {
    for (int d = 0; d < 3; d++)
      if (std::fabs(nodeList(srcIndex,d)-nodeList(tarIndex,d)) > 1e-4*dxMin) 
        return false;
    return true;
  }
}
