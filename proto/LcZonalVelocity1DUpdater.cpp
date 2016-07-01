/**
 * @file	LcZonalVelocity1DUpdater.cpp
 *
 * @brief	
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcZonalVelocity1DUpdater.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  const char *ZonalVelocity1DUpdater::id = "ZonalVelocity1DUpdater";

  ZonalVelocity1DUpdater::ZonalVelocity1DUpdater()
  {
  }

  ZonalVelocity1DUpdater::~ZonalVelocity1DUpdater()
  {
  }

  void
  ZonalVelocity1DUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis");
    else
      throw Lucee::Except("ZonalVelocity1DUpdater::readInput: Must specify element to use using 'basis'");
 
    // should only increments be computed?
    onlyIncrement = false ;
    if (tbl.hasBool("onlyIncrement"))
      onlyIncrement = tbl.getBool("onlyIncrement");
  }

  void
  ZonalVelocity1DUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();
    
// get hold of grid
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();
    // local region to update
    Lucee::Region<1, int> localRgn = grid.getLocalRegion();

    unsigned nlocal = nodalBasis->getNumNodes();

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

    massMatrixInv = Eigen::MatrixXd(nlocal, nlocal);
    massMatrixInv = massMatrix.inverse();

    // Store stiffness matrix
    Lucee::Matrix<double> tempMatrix(nlocal, nlocal);
    nodalBasis->getGradStiffnessMatrix(0, tempMatrix);
    gradStiffnessMatrix = Eigen::MatrixXd(nlocal, nlocal);
    copyLuceeToEigen(tempMatrix, gradStiffnessMatrix);

    // get number of surface quadrature points
    int nSurfQuad = nodalBasis->getNumSurfGaussNodes();

    // temporary variables
    lowerSurfWeights = std::vector<double>(nSurfQuad);
    upperSurfWeights = std::vector<double>(nSurfQuad);

    Lucee::Matrix<double> tempSurfQuad(nSurfQuad, nlocal);
    Lucee::Matrix<double> tempSurfCoords(nSurfQuad, 3);

    // lower surface data
    nodalBasis->getSurfLowerGaussQuadData(0, tempSurfQuad,
      tempSurfCoords, lowerSurfWeights);
    copyLuceeToEigen(tempSurfQuad, lowerSurfInterpMatrix);

    // upper surface data
    nodalBasis->getSurfUpperGaussQuadData(0, tempSurfQuad,
      tempSurfCoords, upperSurfWeights);
    copyLuceeToEigen(tempSurfQuad, upperSurfInterpMatrix);
  }

  Lucee::UpdaterStatus
  ZonalVelocity1DUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();

    const Lucee::Field<1, double>& inpFld = this->getInp<Lucee::Field<1, double> >(0);
    Lucee::Field<1, double>& outFld = this->getOut<Lucee::Field<1, double> >(0);

    double dt = t-this->getCurrTime();
    
    Lucee::ConstFieldPtr<double> inpFldPtr  = inpFld.createConstPtr();
    Lucee::ConstFieldPtr<double> inpFldPtr_r  = inpFld.createConstPtr();
    Lucee::ConstFieldPtr<double> inpFldPtr_l  = inpFld.createConstPtr();
    Lucee::FieldPtr<double> outFldPtr = outFld.createPtr();
    Lucee::FieldPtr<double> outFldPtr_r = outFld.createPtr();
    Lucee::FieldPtr<double> outFldPtr_l = outFld.createPtr();

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

    // Loop over local region cells for volume integral update
    for (int ix = localRgn.getLower(0); ix < localRgn.getUpper(0); ix++)
    {
      inpFld.setPtr(inpFldPtr, ix);
      outFld.setPtr(outFldPtr, ix);
      
      Eigen::VectorXd fAtNodes(nlocal);

      for (int i = 0; i < nlocal; i++)
        fAtNodes(i) = inpFldPtr[i];

      Eigen::VectorXd rhsVector = massMatrixInv*gradStiffnessMatrix*fAtNodes;
      
      // Accumulate updateF to output
      for (int i = 0; i < nlocal; i++)
        outFldPtr[i] = outFldPtr[i] - rhsVector(i);
    }

    // Loop over cell edges
    for (int ix = localRgn.getLower(0); ix < localRgn.getUpper(0) + 1; ix++)
    {
      inpFld.setPtr(inpFldPtr_l, ix-1);
      inpFld.setPtr(inpFldPtr_r, ix);
      outFld.setPtr(outFldPtr_l, ix-1);
      outFld.setPtr(outFldPtr_r, ix);
      
      Eigen::VectorXd fLeftAtNodes(nlocal);
      Eigen::VectorXd fRightAtNodes(nlocal);

      for (int i = 0; i < nlocal; i++)
      {
        fLeftAtNodes(i) = inpFldPtr_l[i];
        fRightAtNodes(i) = inpFldPtr_r[i];
      }

      // Interpolate solution to cell edges
      double fLeftAtSurf = fLeftAtNodes(nlocal-1);//upperSurfInterpMatrix*fLeftAtNodes;
      double fRightAtSurf = fRightAtNodes(0);//lowerSurfInterpMatrix*fRightAtNodes;
      // Take a central flux
      double numFlux = 0.5*(fLeftAtSurf + fRightAtSurf);

      Eigen::VectorXd resultVectorLeft = Eigen::VectorXd::Zero(nlocal);
      Eigen::VectorXd resultVectorRight = Eigen::VectorXd::Zero(nlocal);
      
      resultVectorLeft(nlocal-1) = numFlux;
      resultVectorRight(0) = numFlux;

      Eigen::VectorXd updateVectorLeft = massMatrixInv*resultVectorLeft;
      Eigen::VectorXd updateVectorRight = massMatrixInv*resultVectorRight;

      // Accumulate updateF to output
      for (int i = 0; i < nlocal; i++)
      {
        outFldPtr_l[i] = outFldPtr_l[i] + updateVectorLeft(i);
        outFldPtr_r[i] = outFldPtr_r[i] - updateVectorRight(i);
      }
    }

    return Lucee::UpdaterStatus(true, dt*cfl/cfla);
  }

  void
  ZonalVelocity1DUpdater::declareTypes()
  {
    // Input field
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // Output: some derivative of input field
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }

  void
  ZonalVelocity1DUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
