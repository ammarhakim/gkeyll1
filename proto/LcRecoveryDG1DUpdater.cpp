/**
 * @file	LcRecoveryDG1DUpdater.cpp
 *
 * @brief	Updater to perform recovery DG calculation for 1d problems.
 * Currently supports calculation of second and third derivatives of input field
 * Note that first derivative calculation doesn't work
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcRecoveryDG1DUpdater.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  const char *RecoveryDG1DUpdater::id = "RecoveryDG1DUpdater";

  RecoveryDG1DUpdater::RecoveryDG1DUpdater()
  {
  }

  RecoveryDG1DUpdater::~RecoveryDG1DUpdater()
  {
  }

  void
  RecoveryDG1DUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis");
    else
      throw Lucee::Except("RecoveryDG1DUpdater::readInput: Must specify element to use using 'basis'");
 
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

    if (tbl.hasNumber("polyOrder"))
      polyOrder = tbl.getNumber("polyOrder");
    else
      throw Lucee::Except("RecoveryDG1DUpdater::readInput: Must specify basis function order using 'polyOrder'");

    if (tbl.hasNumber("derivOrder"))
    {
      derivOrder = tbl.getNumber("derivOrder");
    }
    else
      throw Lucee::Except("RecoveryDG1DUpdater::readInput: Must specify derivative order using 'derivOrder'");

    if (polyOrder > 2 || derivOrder > 3)
     throw Lucee::Except("RecoveryDG1DUpdater::readInput: Only polyOrder 1-2 and derivOrder 1-3 supported.");
  }

  void
  RecoveryDG1DUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();
    
// get hold of grid
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();
    // local region to update
    Lucee::Region<1, int> localRgn = grid.getLocalRegion();

    unsigned nlocal = nodalBasis->getNumNodes();
    lowerMat = Eigen::MatrixXd(nlocal, nlocal);
    upperMat = Eigen::MatrixXd(nlocal, nlocal);

    Lucee::RowMajorSequencer<1> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[1];
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
    selfCenter = Eigen::MatrixXd(nlocal, nlocal);
    lowerCenter = Eigen::MatrixXd(nlocal, nlocal);

    double dq[1];
    dq[0] = grid.getDx(0);

    if (polyOrder == 1)
    {
      if (derivOrder == 1)
      {
        selfCenter(0,0)=1/2;
        selfCenter(0,1)=1/2;
        selfCenter(1,0)=-1/2;
        selfCenter(1,1)=-1/2;
        upperCenter(0,0)=0;
        upperCenter(0,1)=0;
        upperCenter(1,0)=1/12;
        upperCenter(1,1)=5/12;
        lowerCenter(0,0)=-5/12;
        lowerCenter(0,1)=-1/12;
        lowerCenter(1,0)=0;
        lowerCenter(1,1)=0;
        upperMat(0,0)=0;
        upperMat(0,1)=0;
        upperMat(1,0)=5/12;
        upperMat(1,1)=1/12;
        lowerMat(0,0)=-1/12;
        lowerMat(0,1)=-5/12;
        lowerMat(1,0)=0;
        lowerMat(1,1)=0;
      }
      else if (derivOrder == 2)
      {
        selfCenter(0,0)=0;
        selfCenter(0,1)=0;
        selfCenter(1,0)=0;
        selfCenter(1,1)=0;
        lowerCenter(0,0)=-13/(6*dq[0]);
        lowerCenter(0,1)=-7/(12*dq[0]);
        lowerCenter(1,0)=5/(12*dq[0]);
        lowerCenter(1,1)=1/(12*dq[0]);
        upperCenter(0,0)=1/(12*dq[0]);
        upperCenter(0,1)=5/(12*dq[0]);
        upperCenter(1,0)=-7/(12*dq[0]);
        upperCenter(1,1)=-13/(6*dq[0]);
        lowerMat(0,0)=5/(12*dq[0]);
        lowerMat(0,1)=4/(3*dq[0]);
        lowerMat(1,0)=1/(12*dq[0]);
        lowerMat(1,1)=5/(12*dq[0]);
        upperMat(0,0)=5/(12*dq[0]);
        upperMat(0,1)=1/(12*dq[0]);
        upperMat(1,0)=4/(3*dq[0]);
        upperMat(1,1)=5/(12*dq[0]);
      }
      else if (derivOrder == 3)
      {
        selfCenter(0,0)=0;
        selfCenter(0,1)=0;
        selfCenter(1,0)=0;
        selfCenter(1,1)=0;
        upperCenter(0,0)=-1/(2*dq[0]*dq[0]);
        upperCenter(0,1)=-7/(4*dq[0]*dq[0]);
        upperCenter(1,0)=3/(2*dq[0]*dq[0]);
        upperCenter(1,1)=3/(4*dq[0]*dq[0]);
        lowerCenter(0,0)=-3/(4*dq[0]*dq[0]);
        lowerCenter(0,1)=-3/(2*dq[0]*dq[0]);
        lowerCenter(1,0)=7/(4*dq[0]*dq[0]);
        lowerCenter(1,1)=1/(2*dq[0]*dq[0]);
        upperMat(0,0)=7/(4*dq[0]*dq[0]);
        upperMat(0,1)=1/(2*dq[0]*dq[0]);
        upperMat(1,0)=-11/(4*dq[0]*dq[0]);
        upperMat(1,1)=1/(2*dq[0]*dq[0]);
        lowerMat(0,0)=-1/(2*dq[0]*dq[0]);
        lowerMat(0,1)=11/(4*dq[0]*dq[0]);
        lowerMat(1,0)=-1/(2*dq[0]*dq[0]);
        lowerMat(1,1)=-7/(4*dq[0]*dq[0]);
      }
    }
    else if (polyOrder == 2)
    {
      if (derivOrder == 1)
      {
        selfCenter(0,0)=1/2;
        selfCenter(0,1)=2/3;
        selfCenter(0,2)=-1/6;
        selfCenter(1,0)=-2/3;
        selfCenter(1,1)=0;
        selfCenter(1,2)=2/3;
        selfCenter(2,0)=1/6;
        selfCenter(2,1)=-2/3;
        selfCenter(2,2)=-1/2;
        upperCenter(0,0)=0;
        upperCenter(0,1)=0;
        upperCenter(0,2)=0;
        upperCenter(1,0)=0;
        upperCenter(1,1)=0;
        upperCenter(1,2)=0;
        upperCenter(2,0)=-3/64;
        upperCenter(2,1)=3/16;
        upperCenter(2,2)=23/64;
        lowerCenter(0,0)=-23/64;
        lowerCenter(0,1)=-3/16;
        lowerCenter(0,2)=3/64;
        lowerCenter(1,0)=0;
        lowerCenter(1,1)=0;
        lowerCenter(1,2)=0;
        lowerCenter(2,0)=0;
        lowerCenter(2,1)=0;
        lowerCenter(2,2)=0;
        upperMat(0,0)=0;
        upperMat(0,1)=0;
        upperMat(0,2)=0;
        upperMat(1,0)=0;
        upperMat(1,1)=0;
        upperMat(1,2)=0;
        upperMat(2,0)=23/64;
        upperMat(2,1)=3/16;
        upperMat(2,2)=-3/64;
        lowerMat(0,0)=3/64;
        lowerMat(0,1)=-3/16;
        lowerMat(0,2)=-23/64;
        lowerMat(1,0)=0;
        lowerMat(1,1)=0;
        lowerMat(1,2)=0;
        lowerMat(2,0)=0;
        lowerMat(2,1)=0;
        lowerMat(2,2)=0;
      }
      else if (derivOrder == 2)
      {
        selfCenter(0,0)=2/(3*dq[0]);
        selfCenter(0,1)=8/(3*dq[0]);
        selfCenter(0,2)=2/(3*dq[0]);
        selfCenter(1,0)=-4/(3*dq[0]);
        selfCenter(1,1)=-16/(3*dq[0]);
        selfCenter(1,2)=-4/(3*dq[0]);
        selfCenter(2,0)=2/(3*dq[0]);
        selfCenter(2,1)=8/(3*dq[0]);
        selfCenter(2,2)=2/(3*dq[0]);
        upperCenter(0,0)=3/(64*dq[0]);
        upperCenter(0,1)=-3/(16*dq[0]);
        upperCenter(0,2)=-23/(64*dq[0]);
        upperCenter(1,0)=-3/(16*dq[0]);
        upperCenter(1,1)=3/(4*dq[0]);
        upperCenter(1,2)=23/(16*dq[0]);
        upperCenter(2,0)=157/(320*dq[0]);
        upperCenter(2,1)=-181/(80*dq[0]);
        upperCenter(2,2)=-1113/(320*dq[0]);
        lowerCenter(0,0)=-1113/(320*dq[0]);
        lowerCenter(0,1)=-181/(80*dq[0]);
        lowerCenter(0,2)=157/(320*dq[0]);
        lowerCenter(1,0)=23/(16*dq[0]);
        lowerCenter(1,1)=3/(4*dq[0]);
        lowerCenter(1,2)=-3/(16*dq[0]);
        lowerCenter(2,0)=-23/(64*dq[0]);
        lowerCenter(2,1)=-3/(16*dq[0]);
        lowerCenter(2,2)=3/(64*dq[0]);
        upperMat(0,0)=-23/(64*dq[0]);
        upperMat(0,1)=-3/(16*dq[0]);
        upperMat(0,2)=3/(64*dq[0]);
        upperMat(1,0)=23/(16*dq[0]);
        upperMat(1,1)=3/(4*dq[0]);
        upperMat(1,2)=-3/(16*dq[0]);
        upperMat(2,0)=423/(320*dq[0]);
        upperMat(2,1)=91/(80*dq[0]);
        upperMat(2,2)=-67/(320*dq[0]);
        lowerMat(0,0)=-67/(320*dq[0]);
        lowerMat(0,1)=91/(80*dq[0]);
        lowerMat(0,2)=423/(320*dq[0]);
        lowerMat(1,0)=-3/(16*dq[0]);
        lowerMat(1,1)=3/(4*dq[0]);
        lowerMat(1,2)=23/(16*dq[0]);
        lowerMat(2,0)=3/(64*dq[0]);
        lowerMat(2,1)=-3/(16*dq[0]);
        lowerMat(2,2)=-23/(64*dq[0]);
      }
      else if (derivOrder == 3)
      {
        selfCenter(0,0)=0;
        selfCenter(0,1)=0;
        selfCenter(0,2)=0;
        selfCenter(1,0)=0;
        selfCenter(1,1)=0;
        selfCenter(1,2)=0;
        selfCenter(2,0)=0;
        selfCenter(2,1)=0;
        selfCenter(2,2)=0;
        upperCenter(0,0)=-43/(80*dq[0]*dq[0]);
        upperCenter(0,1)=49/(20*dq[0]*dq[0]);
        upperCenter(0,2)=307/(80*dq[0]*dq[0]);
        upperCenter(1,0)=71/(40*dq[0]*dq[0]);
        upperCenter(1,1)=-83/(10*dq[0]*dq[0]);
        upperCenter(1,2)=-499/(40*dq[0]*dq[0]);
        upperCenter(2,0)=-89/(80*dq[0]*dq[0]);
        upperCenter(2,1)=187/(20*dq[0]*dq[0]);
        upperCenter(2,2)=401/(80*dq[0]*dq[0]);
        lowerCenter(0,0)=-401/(80*dq[0]*dq[0]);
        lowerCenter(0,1)=-187/(20*dq[0]*dq[0]);
        lowerCenter(0,2)=89/(80*dq[0]*dq[0]);
        lowerCenter(1,0)=499/(40*dq[0]*dq[0]);
        lowerCenter(1,1)=83/(10*dq[0]*dq[0]);
        lowerCenter(1,2)=-71/(40*dq[0]*dq[0]);
        lowerCenter(2,0)=-307/(80*dq[0]*dq[0]);
        lowerCenter(2,1)=-49/(20*dq[0]*dq[0]);
        lowerCenter(2,2)=43/(80*dq[0]*dq[0]);
        upperMat(0,0)=-77/(80*dq[0]*dq[0]);
        upperMat(0,1)=-19/(20*dq[0]*dq[0]);
        upperMat(0,2)=13/(80*dq[0]*dq[0]);
        upperMat(1,0)=269/(40*dq[0]*dq[0]);
        upperMat(1,1)=53/(10*dq[0]*dq[0]);
        upperMat(1,2)=-41/(40*dq[0]*dq[0]);
        upperMat(2,0)=-751/(80*dq[0]*dq[0]);
        upperMat(2,1)=-17/(20*dq[0]*dq[0]);
        upperMat(2,2)=79/(80*dq[0]*dq[0]);
        lowerMat(0,0)=-79/(80*dq[0]*dq[0]);
        lowerMat(0,1)=17/(20*dq[0]*dq[0]);
        lowerMat(0,2)=751/(80*dq[0]*dq[0]);
        lowerMat(1,0)=41/(40*dq[0]*dq[0]);
        lowerMat(1,1)=-53/(10*dq[0]*dq[0]);
        lowerMat(1,2)=-269/(40*dq[0]*dq[0]);
        lowerMat(2,0)=-13/(80*dq[0]*dq[0]);
        lowerMat(2,1)=19/(20*dq[0]*dq[0]);
        lowerMat(2,2)=77/(80*dq[0]*dq[0]);
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
  RecoveryDG1DUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();

    const Lucee::Field<1, double>& inpFld = this->getInp<Lucee::Field<1, double> >(0);
    Lucee::Field<1, double>& outFld = this->getOut<Lucee::Field<1, double> >(0);

    double dt = t-this->getCurrTime();
    
    Lucee::ConstFieldPtr<double> inpFldPtr  = inpFld.createConstPtr();
    Lucee::FieldPtr<double> outFldPtr = outFld.createPtr();

    // check time-step
    double cflm = 1.1*cfl;
    double cfla = 0.0;
    
    outFld = 0.0;

    // local region to index
    Lucee::Region<1, int> localRgn = grid.getLocalRegion();
    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion();

    Lucee::RowMajorSequencer<1> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[1];
    double xc[1];
    seq.fillWithIndex(idx);
    nodalBasis->setIndex(idx);
    unsigned nlocal = nodalBasis->getNumNodes(); 

    // Loop over local region cells
    for (int ix = localRgn.getLower(0); ix < localRgn.getUpper(0); ix++)
    {
      inpFld.setPtr(inpFldPtr, ix);
      outFld.setPtr(outFldPtr, ix);
      
      // Keep track of maximum cfla (v-parallel diffusion)
      cfla = std::max(cfla, dt/(grid.getDx(0)*grid.getDx(0)) );
      if (checkTimeStepSize == true && cfla > cflm)
        return Lucee::UpdaterStatus(false, dt*cfl/cfla);
      
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

      // add in contribution from cells attached to lower/upper faces in vPara
      //if (ix > globalRgn.getLower(0))
      {
        updateF = updateF + lowerCenter*fAtNodes;
        inpFld.setPtr(inpFldPtr, ix-1); // cell attached to lower face
        for (int i = 0; i < nlocal; i++)
          fLowerAtNodes(i) = inpFldPtr[i];
        updateF = updateF + lowerMat*fLowerAtNodes;
      }

      //if (ix < globalRgn.getUpper(0)-1)
      {
        updateF = updateF + upperCenter*fAtNodes;
        inpFld.setPtr(inpFldPtr, ix+1); // cell attached to upper face
        for (int i = 0; i < nlocal; i++)
          fUpperAtNodes(i) = inpFldPtr[i];
        updateF = updateF + upperMat*fUpperAtNodes;
      }

      // Accumulate updateF to output
      for (int i = 0; i < nlocal; i++)
        outFldPtr[i] = outFldPtr[i] + updateF(i);
    }

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
  RecoveryDG1DUpdater::declareTypes()
  {
    // Input field
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // Output: some derivative of input field
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }

  void
  RecoveryDG1DUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
