/**
 * @file	LcRecoveryDG3DUpdater.cpp
 *
 * @brief	Updater to perform recovery DG calculation for 3d problems.
 * Currently supports calculation of second and fourth derivatives of polyOrder = 1 input field
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcRecoveryDG3DUpdater.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  const char *RecoveryDG3DUpdater::id = "RecoveryDG3DUpdater";

  RecoveryDG3DUpdater::RecoveryDG3DUpdater()
  {
  }

  RecoveryDG3DUpdater::~RecoveryDG3DUpdater()
  {
  }

  void
  RecoveryDG3DUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis");
    else
      throw Lucee::Except("RecoveryDG3DUpdater::readInput: Must specify element to use using 'basis'");
 
    // should only increments be computed?
    onlyIncrement = false ;
    if (tbl.hasBool("onlyIncrement"))
      onlyIncrement = tbl.getBool("onlyIncrement");

    // if false, we will ignore cfl checks
    checkTimeStepSize = true;
    // CFL number to control time-step
    if (tbl.hasNumber("cfl"))
      cfl = tbl.getNumber("cfl"); // CFL number
    else
      checkTimeStepSize = false;

    // Factor in front of derivative
    alpha = 1.0;
    if (tbl.hasNumber("alpha"))
      alpha = tbl.getNumber("alpha");

    if (tbl.hasNumber("polyOrder"))
      polyOrder = tbl.getNumber("polyOrder");
    else
      throw Lucee::Except("RecoveryDG3DUpdater::readInput: Must specify basis function order using 'polyOrder'");

    if (tbl.hasNumber("derivOrder"))
    {
      derivOrder = tbl.getNumber("derivOrder");
    }
    else
      throw Lucee::Except("RecoveryDG3DUpdater::readInput: Must specify derivative order using 'derivOrder'");

    if (polyOrder != 1 || (derivOrder != 2 && derivOrder != 4) )
     throw Lucee::Except("RecoveryDG3DUpdater::readInput: Only polyOrder 1 and derivOrder 2 and 4 supported.");
  }

  void
  RecoveryDG3DUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();
    
// get hold of grid
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();
    // local region to update
    Lucee::Region<3, int> localRgn = grid.getLocalRegion();

    unsigned nlocal = nodalBasis->getNumNodes();
    lowerMat = Eigen::MatrixXd(nlocal, nlocal);
    upperMat = Eigen::MatrixXd(nlocal, nlocal);

    Lucee::RowMajorSequencer<3> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[3];
    seq.fillWithIndex(idx);
    nodalBasis->setIndex(idx);

    // pre-multiply each of the matrices by inverse matrix
    Lucee::Matrix<double> massMatrixLucee(nlocal, nlocal);
    nodalBasis->getMassMatrix(massMatrixLucee);
    Eigen::MatrixXd massMatrix(nlocal, nlocal);
    copyLuceeToEigen(massMatrixLucee, massMatrix);
    Eigen::MatrixXd massMatrixInv = massMatrix.inverse();

    // Additional matrices for testing
    upperCenter = Eigen::MatrixXd(nlocal, nlocal);
    selfCenter = Eigen::MatrixXd::Zero(nlocal, nlocal);
    lowerCenter = Eigen::MatrixXd(nlocal, nlocal);

    double dq[3];
    dq[0] = grid.getDx(0);
    dq[1] = grid.getDx(1);
    dq[2] = grid.getDx(2);

    if (polyOrder == 1)
    {
      if (derivOrder == 2)
      {
        upperCenter(0,0)=(dq[1]*dq[2])/(108*dq[0]);
        upperCenter(0,1)=(5*dq[1]*dq[2])/(108*dq[0]);
        upperCenter(0,2)=(dq[1]*dq[2])/(216*dq[0]);
        upperCenter(0,3)=(5*dq[1]*dq[2])/(216*dq[0]);
        upperCenter(0,4)=(dq[1]*dq[2])/(216*dq[0]);
        upperCenter(0,5)=(5*dq[1]*dq[2])/(216*dq[0]);
        upperCenter(0,6)=(dq[1]*dq[2])/(432*dq[0]);
        upperCenter(0,7)=(5*dq[1]*dq[2])/(432*dq[0]);
        upperCenter(1,0)=(-7*dq[1]*dq[2])/(108*dq[0]);
        upperCenter(1,1)=(-13*dq[1]*dq[2])/(54*dq[0]);
        upperCenter(1,2)=(-7*dq[1]*dq[2])/(216*dq[0]);
        upperCenter(1,3)=(-13*dq[1]*dq[2])/(108*dq[0]);
        upperCenter(1,4)=(-7*dq[1]*dq[2])/(216*dq[0]);
        upperCenter(1,5)=(-13*dq[1]*dq[2])/(108*dq[0]);
        upperCenter(1,6)=(-7*dq[1]*dq[2])/(432*dq[0]);
        upperCenter(1,7)=(-13*dq[1]*dq[2])/(216*dq[0]);
        upperCenter(2,0)=(dq[1]*dq[2])/(216*dq[0]);
        upperCenter(2,1)=(5*dq[1]*dq[2])/(216*dq[0]);
        upperCenter(2,2)=(dq[1]*dq[2])/(108*dq[0]);
        upperCenter(2,3)=(5*dq[1]*dq[2])/(108*dq[0]);
        upperCenter(2,4)=(dq[1]*dq[2])/(432*dq[0]);
        upperCenter(2,5)=(5*dq[1]*dq[2])/(432*dq[0]);
        upperCenter(2,6)=(dq[1]*dq[2])/(216*dq[0]);
        upperCenter(2,7)=(5*dq[1]*dq[2])/(216*dq[0]);
        upperCenter(3,0)=(-7*dq[1]*dq[2])/(216*dq[0]);
        upperCenter(3,1)=(-13*dq[1]*dq[2])/(108*dq[0]);
        upperCenter(3,2)=(-7*dq[1]*dq[2])/(108*dq[0]);
        upperCenter(3,3)=(-13*dq[1]*dq[2])/(54*dq[0]);
        upperCenter(3,4)=(-7*dq[1]*dq[2])/(432*dq[0]);
        upperCenter(3,5)=(-13*dq[1]*dq[2])/(216*dq[0]);
        upperCenter(3,6)=(-7*dq[1]*dq[2])/(216*dq[0]);
        upperCenter(3,7)=(-13*dq[1]*dq[2])/(108*dq[0]);
        upperCenter(4,0)=(dq[1]*dq[2])/(216*dq[0]);
        upperCenter(4,1)=(5*dq[1]*dq[2])/(216*dq[0]);
        upperCenter(4,2)=(dq[1]*dq[2])/(432*dq[0]);
        upperCenter(4,3)=(5*dq[1]*dq[2])/(432*dq[0]);
        upperCenter(4,4)=(dq[1]*dq[2])/(108*dq[0]);
        upperCenter(4,5)=(5*dq[1]*dq[2])/(108*dq[0]);
        upperCenter(4,6)=(dq[1]*dq[2])/(216*dq[0]);
        upperCenter(4,7)=(5*dq[1]*dq[2])/(216*dq[0]);
        upperCenter(5,0)=(-7*dq[1]*dq[2])/(216*dq[0]);
        upperCenter(5,1)=(-13*dq[1]*dq[2])/(108*dq[0]);
        upperCenter(5,2)=(-7*dq[1]*dq[2])/(432*dq[0]);
        upperCenter(5,3)=(-13*dq[1]*dq[2])/(216*dq[0]);
        upperCenter(5,4)=(-7*dq[1]*dq[2])/(108*dq[0]);
        upperCenter(5,5)=(-13*dq[1]*dq[2])/(54*dq[0]);
        upperCenter(5,6)=(-7*dq[1]*dq[2])/(216*dq[0]);
        upperCenter(5,7)=(-13*dq[1]*dq[2])/(108*dq[0]);
        upperCenter(6,0)=(dq[1]*dq[2])/(432*dq[0]);
        upperCenter(6,1)=(5*dq[1]*dq[2])/(432*dq[0]);
        upperCenter(6,2)=(dq[1]*dq[2])/(216*dq[0]);
        upperCenter(6,3)=(5*dq[1]*dq[2])/(216*dq[0]);
        upperCenter(6,4)=(dq[1]*dq[2])/(216*dq[0]);
        upperCenter(6,5)=(5*dq[1]*dq[2])/(216*dq[0]);
        upperCenter(6,6)=(dq[1]*dq[2])/(108*dq[0]);
        upperCenter(6,7)=(5*dq[1]*dq[2])/(108*dq[0]);
        upperCenter(7,0)=(-7*dq[1]*dq[2])/(432*dq[0]);
        upperCenter(7,1)=(-13*dq[1]*dq[2])/(216*dq[0]);
        upperCenter(7,2)=(-7*dq[1]*dq[2])/(216*dq[0]);
        upperCenter(7,3)=(-13*dq[1]*dq[2])/(108*dq[0]);
        upperCenter(7,4)=(-7*dq[1]*dq[2])/(216*dq[0]);
        upperCenter(7,5)=(-13*dq[1]*dq[2])/(108*dq[0]);
        upperCenter(7,6)=(-7*dq[1]*dq[2])/(108*dq[0]);
        upperCenter(7,7)=(-13*dq[1]*dq[2])/(54*dq[0]);
        lowerCenter(0,0)=(-13*dq[1]*dq[2])/(54*dq[0]);
        lowerCenter(0,1)=(-7*dq[1]*dq[2])/(108*dq[0]);
        lowerCenter(0,2)=(-13*dq[1]*dq[2])/(108*dq[0]);
        lowerCenter(0,3)=(-7*dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(0,4)=(-13*dq[1]*dq[2])/(108*dq[0]);
        lowerCenter(0,5)=(-7*dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(0,6)=(-13*dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(0,7)=(-7*dq[1]*dq[2])/(432*dq[0]);
        lowerCenter(1,0)=(5*dq[1]*dq[2])/(108*dq[0]);
        lowerCenter(1,1)=(dq[1]*dq[2])/(108*dq[0]);
        lowerCenter(1,2)=(5*dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(1,3)=(dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(1,4)=(5*dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(1,5)=(dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(1,6)=(5*dq[1]*dq[2])/(432*dq[0]);
        lowerCenter(1,7)=(dq[1]*dq[2])/(432*dq[0]);
        lowerCenter(2,0)=(-13*dq[1]*dq[2])/(108*dq[0]);
        lowerCenter(2,1)=(-7*dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(2,2)=(-13*dq[1]*dq[2])/(54*dq[0]);
        lowerCenter(2,3)=(-7*dq[1]*dq[2])/(108*dq[0]);
        lowerCenter(2,4)=(-13*dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(2,5)=(-7*dq[1]*dq[2])/(432*dq[0]);
        lowerCenter(2,6)=(-13*dq[1]*dq[2])/(108*dq[0]);
        lowerCenter(2,7)=(-7*dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(3,0)=(5*dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(3,1)=(dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(3,2)=(5*dq[1]*dq[2])/(108*dq[0]);
        lowerCenter(3,3)=(dq[1]*dq[2])/(108*dq[0]);
        lowerCenter(3,4)=(5*dq[1]*dq[2])/(432*dq[0]);
        lowerCenter(3,5)=(dq[1]*dq[2])/(432*dq[0]);
        lowerCenter(3,6)=(5*dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(3,7)=(dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(4,0)=(-13*dq[1]*dq[2])/(108*dq[0]);
        lowerCenter(4,1)=(-7*dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(4,2)=(-13*dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(4,3)=(-7*dq[1]*dq[2])/(432*dq[0]);
        lowerCenter(4,4)=(-13*dq[1]*dq[2])/(54*dq[0]);
        lowerCenter(4,5)=(-7*dq[1]*dq[2])/(108*dq[0]);
        lowerCenter(4,6)=(-13*dq[1]*dq[2])/(108*dq[0]);
        lowerCenter(4,7)=(-7*dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(5,0)=(5*dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(5,1)=(dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(5,2)=(5*dq[1]*dq[2])/(432*dq[0]);
        lowerCenter(5,3)=(dq[1]*dq[2])/(432*dq[0]);
        lowerCenter(5,4)=(5*dq[1]*dq[2])/(108*dq[0]);
        lowerCenter(5,5)=(dq[1]*dq[2])/(108*dq[0]);
        lowerCenter(5,6)=(5*dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(5,7)=(dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(6,0)=(-13*dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(6,1)=(-7*dq[1]*dq[2])/(432*dq[0]);
        lowerCenter(6,2)=(-13*dq[1]*dq[2])/(108*dq[0]);
        lowerCenter(6,3)=(-7*dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(6,4)=(-13*dq[1]*dq[2])/(108*dq[0]);
        lowerCenter(6,5)=(-7*dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(6,6)=(-13*dq[1]*dq[2])/(54*dq[0]);
        lowerCenter(6,7)=(-7*dq[1]*dq[2])/(108*dq[0]);
        lowerCenter(7,0)=(5*dq[1]*dq[2])/(432*dq[0]);
        lowerCenter(7,1)=(dq[1]*dq[2])/(432*dq[0]);
        lowerCenter(7,2)=(5*dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(7,3)=(dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(7,4)=(5*dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(7,5)=(dq[1]*dq[2])/(216*dq[0]);
        lowerCenter(7,6)=(5*dq[1]*dq[2])/(108*dq[0]);
        lowerCenter(7,7)=(dq[1]*dq[2])/(108*dq[0]);
        upperMat(0,0)=(5*dq[1]*dq[2])/(108*dq[0]);
        upperMat(0,1)=(dq[1]*dq[2])/(108*dq[0]);
        upperMat(0,2)=(5*dq[1]*dq[2])/(216*dq[0]);
        upperMat(0,3)=(dq[1]*dq[2])/(216*dq[0]);
        upperMat(0,4)=(5*dq[1]*dq[2])/(216*dq[0]);
        upperMat(0,5)=(dq[1]*dq[2])/(216*dq[0]);
        upperMat(0,6)=(5*dq[1]*dq[2])/(432*dq[0]);
        upperMat(0,7)=(dq[1]*dq[2])/(432*dq[0]);
        upperMat(1,0)=(4*dq[1]*dq[2])/(27*dq[0]);
        upperMat(1,1)=(5*dq[1]*dq[2])/(108*dq[0]);
        upperMat(1,2)=(2*dq[1]*dq[2])/(27*dq[0]);
        upperMat(1,3)=(5*dq[1]*dq[2])/(216*dq[0]);
        upperMat(1,4)=(2*dq[1]*dq[2])/(27*dq[0]);
        upperMat(1,5)=(5*dq[1]*dq[2])/(216*dq[0]);
        upperMat(1,6)=(dq[1]*dq[2])/(27*dq[0]);
        upperMat(1,7)=(5*dq[1]*dq[2])/(432*dq[0]);
        upperMat(2,0)=(5*dq[1]*dq[2])/(216*dq[0]);
        upperMat(2,1)=(dq[1]*dq[2])/(216*dq[0]);
        upperMat(2,2)=(5*dq[1]*dq[2])/(108*dq[0]);
        upperMat(2,3)=(dq[1]*dq[2])/(108*dq[0]);
        upperMat(2,4)=(5*dq[1]*dq[2])/(432*dq[0]);
        upperMat(2,5)=(dq[1]*dq[2])/(432*dq[0]);
        upperMat(2,6)=(5*dq[1]*dq[2])/(216*dq[0]);
        upperMat(2,7)=(dq[1]*dq[2])/(216*dq[0]);
        upperMat(3,0)=(2*dq[1]*dq[2])/(27*dq[0]);
        upperMat(3,1)=(5*dq[1]*dq[2])/(216*dq[0]);
        upperMat(3,2)=(4*dq[1]*dq[2])/(27*dq[0]);
        upperMat(3,3)=(5*dq[1]*dq[2])/(108*dq[0]);
        upperMat(3,4)=(dq[1]*dq[2])/(27*dq[0]);
        upperMat(3,5)=(5*dq[1]*dq[2])/(432*dq[0]);
        upperMat(3,6)=(2*dq[1]*dq[2])/(27*dq[0]);
        upperMat(3,7)=(5*dq[1]*dq[2])/(216*dq[0]);
        upperMat(4,0)=(5*dq[1]*dq[2])/(216*dq[0]);
        upperMat(4,1)=(dq[1]*dq[2])/(216*dq[0]);
        upperMat(4,2)=(5*dq[1]*dq[2])/(432*dq[0]);
        upperMat(4,3)=(dq[1]*dq[2])/(432*dq[0]);
        upperMat(4,4)=(5*dq[1]*dq[2])/(108*dq[0]);
        upperMat(4,5)=(dq[1]*dq[2])/(108*dq[0]);
        upperMat(4,6)=(5*dq[1]*dq[2])/(216*dq[0]);
        upperMat(4,7)=(dq[1]*dq[2])/(216*dq[0]);
        upperMat(5,0)=(2*dq[1]*dq[2])/(27*dq[0]);
        upperMat(5,1)=(5*dq[1]*dq[2])/(216*dq[0]);
        upperMat(5,2)=(dq[1]*dq[2])/(27*dq[0]);
        upperMat(5,3)=(5*dq[1]*dq[2])/(432*dq[0]);
        upperMat(5,4)=(4*dq[1]*dq[2])/(27*dq[0]);
        upperMat(5,5)=(5*dq[1]*dq[2])/(108*dq[0]);
        upperMat(5,6)=(2*dq[1]*dq[2])/(27*dq[0]);
        upperMat(5,7)=(5*dq[1]*dq[2])/(216*dq[0]);
        upperMat(6,0)=(5*dq[1]*dq[2])/(432*dq[0]);
        upperMat(6,1)=(dq[1]*dq[2])/(432*dq[0]);
        upperMat(6,2)=(5*dq[1]*dq[2])/(216*dq[0]);
        upperMat(6,3)=(dq[1]*dq[2])/(216*dq[0]);
        upperMat(6,4)=(5*dq[1]*dq[2])/(216*dq[0]);
        upperMat(6,5)=(dq[1]*dq[2])/(216*dq[0]);
        upperMat(6,6)=(5*dq[1]*dq[2])/(108*dq[0]);
        upperMat(6,7)=(dq[1]*dq[2])/(108*dq[0]);
        upperMat(7,0)=(dq[1]*dq[2])/(27*dq[0]);
        upperMat(7,1)=(5*dq[1]*dq[2])/(432*dq[0]);
        upperMat(7,2)=(2*dq[1]*dq[2])/(27*dq[0]);
        upperMat(7,3)=(5*dq[1]*dq[2])/(216*dq[0]);
        upperMat(7,4)=(2*dq[1]*dq[2])/(27*dq[0]);
        upperMat(7,5)=(5*dq[1]*dq[2])/(216*dq[0]);
        upperMat(7,6)=(4*dq[1]*dq[2])/(27*dq[0]);
        upperMat(7,7)=(5*dq[1]*dq[2])/(108*dq[0]);
        lowerMat(0,0)=(5*dq[1]*dq[2])/(108*dq[0]);
        lowerMat(0,1)=(4*dq[1]*dq[2])/(27*dq[0]);
        lowerMat(0,2)=(5*dq[1]*dq[2])/(216*dq[0]);
        lowerMat(0,3)=(2*dq[1]*dq[2])/(27*dq[0]);
        lowerMat(0,4)=(5*dq[1]*dq[2])/(216*dq[0]);
        lowerMat(0,5)=(2*dq[1]*dq[2])/(27*dq[0]);
        lowerMat(0,6)=(5*dq[1]*dq[2])/(432*dq[0]);
        lowerMat(0,7)=(dq[1]*dq[2])/(27*dq[0]);
        lowerMat(1,0)=(dq[1]*dq[2])/(108*dq[0]);
        lowerMat(1,1)=(5*dq[1]*dq[2])/(108*dq[0]);
        lowerMat(1,2)=(dq[1]*dq[2])/(216*dq[0]);
        lowerMat(1,3)=(5*dq[1]*dq[2])/(216*dq[0]);
        lowerMat(1,4)=(dq[1]*dq[2])/(216*dq[0]);
        lowerMat(1,5)=(5*dq[1]*dq[2])/(216*dq[0]);
        lowerMat(1,6)=(dq[1]*dq[2])/(432*dq[0]);
        lowerMat(1,7)=(5*dq[1]*dq[2])/(432*dq[0]);
        lowerMat(2,0)=(5*dq[1]*dq[2])/(216*dq[0]);
        lowerMat(2,1)=(2*dq[1]*dq[2])/(27*dq[0]);
        lowerMat(2,2)=(5*dq[1]*dq[2])/(108*dq[0]);
        lowerMat(2,3)=(4*dq[1]*dq[2])/(27*dq[0]);
        lowerMat(2,4)=(5*dq[1]*dq[2])/(432*dq[0]);
        lowerMat(2,5)=(dq[1]*dq[2])/(27*dq[0]);
        lowerMat(2,6)=(5*dq[1]*dq[2])/(216*dq[0]);
        lowerMat(2,7)=(2*dq[1]*dq[2])/(27*dq[0]);
        lowerMat(3,0)=(dq[1]*dq[2])/(216*dq[0]);
        lowerMat(3,1)=(5*dq[1]*dq[2])/(216*dq[0]);
        lowerMat(3,2)=(dq[1]*dq[2])/(108*dq[0]);
        lowerMat(3,3)=(5*dq[1]*dq[2])/(108*dq[0]);
        lowerMat(3,4)=(dq[1]*dq[2])/(432*dq[0]);
        lowerMat(3,5)=(5*dq[1]*dq[2])/(432*dq[0]);
        lowerMat(3,6)=(dq[1]*dq[2])/(216*dq[0]);
        lowerMat(3,7)=(5*dq[1]*dq[2])/(216*dq[0]);
        lowerMat(4,0)=(5*dq[1]*dq[2])/(216*dq[0]);
        lowerMat(4,1)=(2*dq[1]*dq[2])/(27*dq[0]);
        lowerMat(4,2)=(5*dq[1]*dq[2])/(432*dq[0]);
        lowerMat(4,3)=(dq[1]*dq[2])/(27*dq[0]);
        lowerMat(4,4)=(5*dq[1]*dq[2])/(108*dq[0]);
        lowerMat(4,5)=(4*dq[1]*dq[2])/(27*dq[0]);
        lowerMat(4,6)=(5*dq[1]*dq[2])/(216*dq[0]);
        lowerMat(4,7)=(2*dq[1]*dq[2])/(27*dq[0]);
        lowerMat(5,0)=(dq[1]*dq[2])/(216*dq[0]);
        lowerMat(5,1)=(5*dq[1]*dq[2])/(216*dq[0]);
        lowerMat(5,2)=(dq[1]*dq[2])/(432*dq[0]);
        lowerMat(5,3)=(5*dq[1]*dq[2])/(432*dq[0]);
        lowerMat(5,4)=(dq[1]*dq[2])/(108*dq[0]);
        lowerMat(5,5)=(5*dq[1]*dq[2])/(108*dq[0]);
        lowerMat(5,6)=(dq[1]*dq[2])/(216*dq[0]);
        lowerMat(5,7)=(5*dq[1]*dq[2])/(216*dq[0]);
        lowerMat(6,0)=(5*dq[1]*dq[2])/(432*dq[0]);
        lowerMat(6,1)=(dq[1]*dq[2])/(27*dq[0]);
        lowerMat(6,2)=(5*dq[1]*dq[2])/(216*dq[0]);
        lowerMat(6,3)=(2*dq[1]*dq[2])/(27*dq[0]);
        lowerMat(6,4)=(5*dq[1]*dq[2])/(216*dq[0]);
        lowerMat(6,5)=(2*dq[1]*dq[2])/(27*dq[0]);
        lowerMat(6,6)=(5*dq[1]*dq[2])/(108*dq[0]);
        lowerMat(6,7)=(4*dq[1]*dq[2])/(27*dq[0]);
        lowerMat(7,0)=(dq[1]*dq[2])/(432*dq[0]);
        lowerMat(7,1)=(5*dq[1]*dq[2])/(432*dq[0]);
        lowerMat(7,2)=(dq[1]*dq[2])/(216*dq[0]);
        lowerMat(7,3)=(5*dq[1]*dq[2])/(216*dq[0]);
        lowerMat(7,4)=(dq[1]*dq[2])/(216*dq[0]);
        lowerMat(7,5)=(5*dq[1]*dq[2])/(216*dq[0]);
        lowerMat(7,6)=(dq[1]*dq[2])/(108*dq[0]);
        lowerMat(7,7)=(5*dq[1]*dq[2])/(108*dq[0]);
      }
      else if (derivOrder == 4)
      {
        upperCenter(0,0)=(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(0,1)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(0,2)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperCenter(0,3)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperCenter(0,4)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperCenter(0,5)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperCenter(0,6)=(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        upperCenter(0,7)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        upperCenter(1,0)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(1,1)=(16*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(1,2)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperCenter(1,3)=(8*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(1,4)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperCenter(1,5)=(8*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(1,6)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        upperCenter(1,7)=(4*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(2,0)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperCenter(2,1)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperCenter(2,2)=(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(2,3)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(2,4)=(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        upperCenter(2,5)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        upperCenter(2,6)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperCenter(2,7)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperCenter(3,0)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperCenter(3,1)=(8*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(3,2)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(3,3)=(16*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(3,4)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        upperCenter(3,5)=(4*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(3,6)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperCenter(3,7)=(8*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(4,0)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperCenter(4,1)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperCenter(4,2)=(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        upperCenter(4,3)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        upperCenter(4,4)=(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(4,5)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(4,6)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperCenter(4,7)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperCenter(5,0)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperCenter(5,1)=(8*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(5,2)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        upperCenter(5,3)=(4*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(5,4)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(5,5)=(16*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(5,6)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperCenter(5,7)=(8*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(6,0)=(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        upperCenter(6,1)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        upperCenter(6,2)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperCenter(6,3)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperCenter(6,4)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperCenter(6,5)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperCenter(6,6)=(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(6,7)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(7,0)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        upperCenter(7,1)=(4*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(7,2)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperCenter(7,3)=(8*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(7,4)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperCenter(7,5)=(8*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(7,6)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperCenter(7,7)=(16*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(0,0)=(16*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(0,1)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(0,2)=(8*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(0,3)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerCenter(0,4)=(8*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(0,5)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerCenter(0,6)=(4*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(0,7)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        lowerCenter(1,0)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(1,1)=(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(1,2)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerCenter(1,3)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerCenter(1,4)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerCenter(1,5)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerCenter(1,6)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        lowerCenter(1,7)=(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        lowerCenter(2,0)=(8*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(2,1)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerCenter(2,2)=(16*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(2,3)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(2,4)=(4*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(2,5)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        lowerCenter(2,6)=(8*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(2,7)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerCenter(3,0)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerCenter(3,1)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerCenter(3,2)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(3,3)=(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(3,4)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        lowerCenter(3,5)=(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        lowerCenter(3,6)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerCenter(3,7)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerCenter(4,0)=(8*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(4,1)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerCenter(4,2)=(4*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(4,3)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        lowerCenter(4,4)=(16*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(4,5)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(4,6)=(8*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(4,7)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerCenter(5,0)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerCenter(5,1)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerCenter(5,2)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        lowerCenter(5,3)=(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        lowerCenter(5,4)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(5,5)=(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(5,6)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerCenter(5,7)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerCenter(6,0)=(4*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(6,1)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        lowerCenter(6,2)=(8*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(6,3)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerCenter(6,4)=(8*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(6,5)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerCenter(6,6)=(16*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(6,7)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(7,0)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        lowerCenter(7,1)=(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        lowerCenter(7,2)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerCenter(7,3)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerCenter(7,4)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerCenter(7,5)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerCenter(7,6)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerCenter(7,7)=(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperMat(0,0)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperMat(0,1)=(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperMat(0,2)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(0,3)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(0,4)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(0,5)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(0,6)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        upperMat(0,7)=(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        upperMat(1,0)=(-14*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperMat(1,1)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperMat(1,2)=(-7*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperMat(1,3)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(1,4)=(-7*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperMat(1,5)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(1,6)=(-7*dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(1,7)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        upperMat(2,0)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(2,1)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(2,2)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperMat(2,3)=(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperMat(2,4)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        upperMat(2,5)=(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        upperMat(2,6)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(2,7)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(3,0)=(-7*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperMat(3,1)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(3,2)=(-14*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperMat(3,3)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperMat(3,4)=(-7*dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(3,5)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        upperMat(3,6)=(-7*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperMat(3,7)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(4,0)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(4,1)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(4,2)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        upperMat(4,3)=(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        upperMat(4,4)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperMat(4,5)=(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperMat(4,6)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(4,7)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(5,0)=(-7*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperMat(5,1)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(5,2)=(-7*dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(5,3)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        upperMat(5,4)=(-14*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperMat(5,5)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperMat(5,6)=(-7*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperMat(5,7)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(6,0)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        upperMat(6,1)=(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        upperMat(6,2)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(6,3)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(6,4)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(6,5)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(6,6)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperMat(6,7)=(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperMat(7,0)=(-7*dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(7,1)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        upperMat(7,2)=(-7*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperMat(7,3)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(7,4)=(-7*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperMat(7,5)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        upperMat(7,6)=(-14*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        upperMat(7,7)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerMat(0,0)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerMat(0,1)=(-14*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerMat(0,2)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(0,3)=(-7*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerMat(0,4)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(0,5)=(-7*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerMat(0,6)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        lowerMat(0,7)=(-7*dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(1,0)=(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerMat(1,1)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerMat(1,2)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(1,3)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(1,4)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(1,5)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(1,6)=(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        lowerMat(1,7)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        lowerMat(2,0)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(2,1)=(-7*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerMat(2,2)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerMat(2,3)=(-14*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerMat(2,4)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        lowerMat(2,5)=(-7*dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(2,6)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(2,7)=(-7*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerMat(3,0)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(3,1)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(3,2)=(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerMat(3,3)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerMat(3,4)=(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        lowerMat(3,5)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        lowerMat(3,6)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(3,7)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(4,0)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(4,1)=(-7*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerMat(4,2)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        lowerMat(4,3)=(-7*dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(4,4)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerMat(4,5)=(-14*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerMat(4,6)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(4,7)=(-7*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerMat(5,0)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(5,1)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(5,2)=(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        lowerMat(5,3)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        lowerMat(5,4)=(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerMat(5,5)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerMat(5,6)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(5,7)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(6,0)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        lowerMat(6,1)=(-7*dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(6,2)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(6,3)=(-7*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerMat(6,4)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(6,5)=(-7*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerMat(6,6)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerMat(6,7)=(-14*dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerMat(7,0)=(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        lowerMat(7,1)=-(dq[1]*dq[2])/(36*dq[0]*dq[0]*dq[0]);
        lowerMat(7,2)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(7,3)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(7,4)=(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(7,5)=-(dq[1]*dq[2])/(18*dq[0]*dq[0]*dq[0]);
        lowerMat(7,6)=(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
        lowerMat(7,7)=-(dq[1]*dq[2])/(9*dq[0]*dq[0]*dq[0]);
      }
    }

    // Multiply relevant matrices by massMatrixInv
    selfCenter = massMatrixInv*selfCenter;
    lowerCenter = massMatrixInv*lowerCenter;
    upperCenter = massMatrixInv*upperCenter;
    lowerMat = massMatrixInv*lowerMat;
    upperMat = massMatrixInv*upperMat;
  }

  Lucee::UpdaterStatus
  RecoveryDG3DUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    const Lucee::Field<3, double>& inpFld = this->getInp<Lucee::Field<3, double> >(0);
    Lucee::Field<3, double>& outFld = this->getOut<Lucee::Field<3, double> >(0);

    double dt = t-this->getCurrTime();
    
    Lucee::ConstFieldPtr<double> inpFldPtr  = inpFld.createConstPtr();
    Lucee::FieldPtr<double> outFldPtr = outFld.createPtr();

    // check time-step
    double cflm = 1.1*cfl;
    double cfla = 0.0;
    
    outFld = 0.0;

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

    seq.reset();

    while(seq.step())
    {
      seq.fillWithIndex(idx);
      inpFld.setPtr(inpFldPtr, idx);
      outFld.setPtr(outFldPtr, idx);
      
      // Keep track of maximum cfla (y-space diffusion)
      if (derivOrder == 2)
        cfla = std::max(cfla, std::fabs(alpha)*4*dt/(grid.getDx(0)*grid.getDx(0)) );
      else if (derivOrder == 4)
        cfla = std::max(cfla, std::fabs(alpha)*16*dt/(grid.getDx(0)*grid.getDx(0)*grid.getDx(0)*grid.getDx(0)));
      
      //idx[0] = ix;
      //grid.setIndex(idx);
      //grid.getCentroid(xc);

      Eigen::VectorXd fAtNodes(nlocal);

      for (int i = 0; i < nlocal; i++)
        fAtNodes(i) = inpFldPtr[i];

      Eigen::VectorXd updateF = selfCenter*fAtNodes;
      
      // Data from neighboring cells
      Eigen::VectorXd fLowerAtNodes(nlocal);
      Eigen::VectorXd fUpperAtNodes(nlocal);

      // add in contribution from cells attached to lower/upper faces in y
      //if (ix > globalRgn.getLower(0))
      {
        updateF = updateF + lowerCenter*fAtNodes;
        inpFld.setPtr(inpFldPtr, idx[0]-1, idx[1], idx[2]); // cell attached to lower face
        for (int i = 0; i < nlocal; i++)
          fLowerAtNodes(i) = inpFldPtr[i];
        updateF = updateF + lowerMat*fLowerAtNodes;
      }

      //if (ix < globalRgn.getUpper(0)-1)
      {
        updateF = updateF + upperCenter*fAtNodes;
        inpFld.setPtr(inpFldPtr, idx[0]+1, idx[1], idx[2]); // cell attached to upper face
        for (int i = 0; i < nlocal; i++)
          fUpperAtNodes(i) = inpFldPtr[i];
        updateF = updateF + upperMat*fUpperAtNodes;
      }

      // Accumulate updateF to output
      for (int i = 0; i < nlocal; i++)
        outFldPtr[i] = outFldPtr[i] + alpha*updateF(i);
    }

    // Check time step criteria
    if (checkTimeStepSize == true && cfla > cflm)
      return Lucee::UpdaterStatus(false, dt*cfl/cfla);

    seq.reset();
    // Final sweep, update solution with forward Euler step

    while (seq.step())
    {
      seq.fillWithIndex(idx);
      outFld.setPtr(outFldPtr, idx);
      grid.setIndex(idx);

      if (onlyIncrement == false)
      {
        inpFld.setPtr(inpFldPtr, idx);
        for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
          outFldPtr[nodeIndex] = inpFldPtr[nodeIndex] + dt*outFldPtr[nodeIndex];
      }
    }

    return Lucee::UpdaterStatus(true, dt*cfl/cfla);
  }

  void
  RecoveryDG3DUpdater::declareTypes()
  {
    // Input field
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Output: some derivative of input field
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
  }

  void
  RecoveryDG3DUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
}
