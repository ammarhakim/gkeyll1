/**
 * @file	LcSOLWeightedProjectionCalc.cpp
 *
 * @brief	Projects the product of B*g*f onto a 3d field, where g and f are 5d fields and g is a 3d field
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLWeightedProjectionCalc.h>

namespace Lucee
{
  const char *SOLWeightedProjectionCalc::id = "SOLWeightedProjectionCalc";

  SOLWeightedProjectionCalc::SOLWeightedProjectionCalc()
  {
  }

  SOLWeightedProjectionCalc::~SOLWeightedProjectionCalc()
  {
    delete result;
  }

  void
  SOLWeightedProjectionCalc::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except("SOLWeightedProjectionCalc::readInput: Must specify element to use using 'basis5d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("SOLWeightedProjectionCalc::readInput: Must specify element to use using 'basis3d'");

    if (tbl.hasNumber("scaleFactor"))
      scaleFactor = tbl.getNumber("scaleFactor");
    else scaleFactor = 1.0;
  }

  void
  SOLWeightedProjectionCalc::initialize()
  {
    UpdaterIfc::initialize();

    unsigned nlocal5d = nodalBasis5d->getNumNodes();
    unsigned nlocal3d = nodalBasis3d->getNumNodes();

    // get volume interpolation matrices for 3d element
    int nVolQuad3d = nodalBasis3d->getNumGaussNodes();
    std::vector<double> gaussWeights3d(nVolQuad3d);
    Lucee::Matrix<double> tempVolQuad3d(nVolQuad3d, nlocal3d);
    Lucee::Matrix<double> tempVolCoords3d(nVolQuad3d, 3);

    nodalBasis3d->getGaussQuadData(tempVolQuad3d, tempVolCoords3d, gaussWeights3d);

    interpMatrix3d = Eigen::MatrixXd(nVolQuad3d, nlocal3d);
    copyLuceeToEigen(tempVolQuad3d, interpMatrix3d);

    // get volume interpolation matrices for 5d element
    int nVolQuad5d = nodalBasis5d->getNumGaussNodes();
    gaussWeights5d = std::vector<double>(nVolQuad5d);
    Lucee::Matrix<double> tempVolQuad5d(nVolQuad5d, nlocal5d);
    Lucee::Matrix<double> tempVolCoords5d(nVolQuad5d, 5);

    nodalBasis5d->getGaussQuadData(tempVolQuad5d, tempVolCoords5d, gaussWeights5d);

    interpMatrix5d = Eigen::MatrixXd(nVolQuad5d, nlocal5d);
    copyLuceeToEigen(tempVolQuad5d, interpMatrix5d);

    // Get and store inverse of 3d mass matrix
    Lucee::Matrix<double> tempMassMatrix3d(nlocal3d, nlocal3d);
    nodalBasis3d->getMassMatrix(tempMassMatrix3d);
    Eigen::MatrixXd massMatrix3d(nlocal3d, nlocal3d);
    copyLuceeToEigen(tempMassMatrix3d, massMatrix3d);

    massMatrixInv3d = massMatrix3d.inverse();

    // get hold of grid
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // local region to update
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();

    int lower[3], upper[3];
    lower[0] = localRgn.getLower(0); upper[0] = localRgn.getUpper(0);
    lower[1] = localRgn.getLower(1); upper[1] = localRgn.getUpper(1);
    lower[2] = localRgn.getLower(2); upper[2] = localRgn.getUpper(2);
    Lucee::Region<3, int> local3D(lower, upper);
    int lg[3] = {1,1,1}, ug[3] = {1,1,1}; // Assume there is a 1 layer ghost cell
    // allocate space for storing local result calculation
    result = new Lucee::Field<3, double>(local3D, nlocal3d, lg, ug);
  }

  Lucee::UpdaterStatus
  SOLWeightedProjectionCalc::update(double t)
  {
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // Distribution function
    const Lucee::Field<5, double>& distfIn = this->getInp<Lucee::Field<5, double> >(0);
    // Field to multiply distribution function by before integrating
    const Lucee::Field<5, double>& distgIn = this->getInp<Lucee::Field<5, double> >(1);
    // Magnetic field in 3d
    const Lucee::Field<3, double>& bFieldIn = this->getInp<Lucee::Field<3, double> >(2);
    // Output field
    Lucee::Field<3, double>& integratedField = this->getOut<Lucee::Field<3, double> >(0);
    
    // create duplicate to store local moments
    //Lucee::Field<3, double> result = momentOut.duplicate();
    // clear out contents of output field
    (*result) = 0.0;

    integratedField= 0.0; // clear out current contents

    Lucee::Region<5, int> globalRgn = grid.getGlobalRegion();
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();
    
    Lucee::ConstFieldPtr<double> distfInPtr = distfIn.createConstPtr();
    Lucee::ConstFieldPtr<double> distgInPtr = distgIn.createConstPtr();
    Lucee::ConstFieldPtr<double> bFieldInPtr = bFieldIn.createConstPtr();
    Lucee::FieldPtr<double> integratedFieldPtr = integratedField.createPtr(); // Output pointer
    Lucee::FieldPtr<double> resultPtr = result->createPtr();

    unsigned nlocal5d = nodalBasis5d->getNumNodes();
    unsigned nlocal3d = nodalBasis3d->getNumNodes();
    int nVolQuad5d = nodalBasis5d->getNumGaussNodes();
    int nVolQuad3d = nodalBasis3d->getNumGaussNodes();

    int idx[5];

    Lucee::RowMajorSequencer<5> seq(localRgn);
    
    Eigen::VectorXd distfVec(nlocal5d);
    Eigen::VectorXd distgVec(nlocal5d);
    Eigen::VectorXd bFieldVec(nlocal3d);

    Eigen::VectorXd distfAtQuad(nVolQuad5d);
    Eigen::VectorXd distgAtQuad(nVolQuad5d);
    Eigen::VectorXd bFieldAtQuad(nVolQuad3d);
    
    Eigen::VectorXd rhsIntegrals(nlocal3d);
    Eigen::VectorXd solutionVec(nlocal3d);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      //result->setPtr(resultPtr, idx[0], idx[1], idx[2]);
      integratedField.setPtr(integratedFieldPtr, idx[0], idx[1], idx[2]);

      distfIn.setPtr(distfInPtr, idx);
      distgIn.setPtr(distgInPtr, idx);
      bFieldIn.setPtr(bFieldInPtr, idx[0], idx[1], idx[2]);

      for (int nodeIndex = 0; nodeIndex < nlocal5d; nodeIndex++)
      {
        distfVec(nodeIndex) = distfInPtr[nodeIndex];
        distgVec(nodeIndex) = distgInPtr[nodeIndex];
      }

      for (int nodeIndex = 0; nodeIndex < nlocal3d; nodeIndex++)
        bFieldVec(nodeIndex) = bFieldInPtr[nodeIndex];

      // Compute fields at quadrature points
      distfAtQuad = interpMatrix5d*distfVec;
      distgAtQuad = interpMatrix5d*distgVec;
      bFieldAtQuad = interpMatrix3d*bFieldVec;
      
      for (int basisIndex = 0; basisIndex < nlocal3d; basisIndex++)
      {
        double integrationResult = 0.0;
        // Compute 5d integration
        for (int quadIndex = 0; quadIndex < nVolQuad5d; quadIndex++)
        {
          integrationResult += gaussWeights5d[quadIndex]*interpMatrix3d(quadIndex % nVolQuad3d, basisIndex)*
            scaleFactor*distfAtQuad(quadIndex)*distgAtQuad(quadIndex)*bFieldAtQuad(quadIndex % nVolQuad3d);
        }
        rhsIntegrals(basisIndex) = integrationResult;

        /*if (std::isnan(integrationResult))
        {
          std::cout << "SOLWeightedProjectionCalc: Integration result nan" << std::endl;
          for (int nodeIndex = 0; nodeIndex < nlocal5d; nodeIndex++)
          {
            std::cout << "distfInPtr[" << nodeIndex << "] = " << distfInPtr[nodeIndex] << std::endl;
          }
        }*/
      }
      // Calculate solution
      solutionVec = massMatrixInv3d*rhsIntegrals;
      
      //for (int nodeIndex = 0; nodeIndex < nlocal3d; nodeIndex++)
      //  resultPtr[nodeIndex] += solutionVec(nodeIndex);

      for (int nodeIndex = 0; nodeIndex < nlocal3d; nodeIndex++)
        integratedFieldPtr[nodeIndex] += solutionVec(nodeIndex);
    }

    /*
    int localPositionCells = localRgn.getShape(0)*localRgn.getShape(1)*localRgn.getShape(2);
    std::vector<double> localResult(localPositionCells*nlocal3d);

    int idx3d[3];
    int lower3d[] = {localRgn.getLower(0), localRgn.getLower(1), localRgn.getLower(2)};
    int upper3d[] = {localRgn.getUpper(0), localRgn.getUpper(1), localRgn.getUpper(2)};
    Lucee::Region<3, int> rgn3d(lower3d, upper3d);
    Lucee::RowMajorIndexer<3> idxr(rgn3d);
    Lucee::RowMajorSequencer<3> seq3d(rgn3d);
    // Loop over each 'local' position space cell
    while(seq3d.step())
    {
      seq3d.fillWithIndex(idx3d);
      int cellIndex = idxr.getIndex(idx3d);

      result->setPtr(resultPtr, idx3d);
      // copy data to vector
      for (int i = 0; i < nlocal3d; i++)
        localResult[cellIndex*nlocal3d+i] = resultPtr[i];
    }

    // Above loop computes results on local phase-space domain. We need to
    // sum across velocity space to get total result on configuration
    // space.
    std::vector<double> reducedMoment(localPositionCells*nlocal3d);
    // we need to get result communicator of field as updater's result
    // communicator is same as its grid's result communicator. In this
    // case, grid is phase-space grid, which is not what we want.
    TxCommBase *momComm = integratedField.getMomComm();
    unsigned xsize = localPositionCells*nlocal3d; // amount to communicate
    momComm->allreduce(xsize, localResult, reducedMoment, TX_SUM);

    seq3d.reset();
    // Copy reducedMoment to output field
    while(seq3d.step())
    {
      seq3d.fillWithIndex(idx3d);
      int cellIndex = idxr.getIndex(idx3d);
      integratedField.setPtr(integratedFieldPtr, idx3d);
      // copy data to vector
      for (int i = 0; i < nlocal3d; i++)
        integratedFieldPtr[i] = reducedMoment[cellIndex*nlocal3d+i];
    }*/
   
    return Lucee::UpdaterStatus();
  }

  void
  SOLWeightedProjectionCalc::declareTypes()
  {
    // Input: distribution function (f)
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Input: field to multiply with the distribution function (g)
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Input: 3d field of magnetic field (B)
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Output: 3d field containing projection of f*g*B
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
  }

  void
  SOLWeightedProjectionCalc::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
