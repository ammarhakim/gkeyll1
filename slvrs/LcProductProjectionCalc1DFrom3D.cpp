/**
 * @file	LcProductProjectionCalc1DFrom3D.cpp
 *
 * @brief	Projects the product of g*f onto a 1d field, where g and f are 3d fields
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcProductProjectionCalc1DFrom3D.h>

namespace Lucee
{
  const char *ProductProjectionCalc1DFrom3D::id = "ProductProjectionCalc1DFrom3D";

  ProductProjectionCalc1DFrom3D::ProductProjectionCalc1DFrom3D()
  {
    resultLocal = 0;
  }

  ProductProjectionCalc1DFrom3D::~ProductProjectionCalc1DFrom3D()
  {
    delete resultLocal;
  }

  void
  ProductProjectionCalc1DFrom3D::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("ProductProjectionCalc1DFrom3D::readInput: Must specify element to use using 'basis3d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis1d"))
      nodalBasis1d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis1d");
    else
      throw Lucee::Except("ProductProjectionCalc1DFrom3D::readInput: Must specify element to use using 'basis1d'");

    if (tbl.hasNumber("scaleFactor"))
      scaleFactor = tbl.getNumber("scaleFactor");
    else scaleFactor = 1.0;
  }

  void
  ProductProjectionCalc1DFrom3D::initialize()
  {
    UpdaterIfc::initialize();

    unsigned nlocal3d = nodalBasis3d->getNumNodes();
    unsigned nlocal1d = nodalBasis1d->getNumNodes();

    // get volume interpolation matrices for 1d element
    int nVolQuad1d = nodalBasis1d->getNumGaussNodes();
    std::vector<double> gaussWeights1d(nVolQuad1d);
    Lucee::Matrix<double> tempVolQuad1d(nVolQuad1d, nlocal1d);
    Lucee::Matrix<double> tempVolCoords1d(nVolQuad1d, 3);

    nodalBasis1d->getGaussQuadData(tempVolQuad1d, tempVolCoords1d, gaussWeights1d);

    interpMatrix1d = Eigen::MatrixXd(nVolQuad1d, nlocal1d);
    copyLuceeToEigen(tempVolQuad1d, interpMatrix1d);

    // get volume interpolation matrices for 3d element
    int nVolQuad3d = nodalBasis3d->getNumGaussNodes();
    gaussWeights3d = std::vector<double>(nVolQuad3d);
    Lucee::Matrix<double> tempVolQuad3d(nVolQuad3d, nlocal3d);
    Lucee::Matrix<double> tempVolCoords3d(nVolQuad3d, 3);

    nodalBasis3d->getGaussQuadData(tempVolQuad3d, tempVolCoords3d, gaussWeights3d);

    interpMatrix3d = Eigen::MatrixXd(nVolQuad3d, nlocal3d);
    copyLuceeToEigen(tempVolQuad3d, interpMatrix3d);

    // Get and store inverse of 1d mass matrix
    Lucee::Matrix<double> tempMassMatrix1d(nlocal1d, nlocal1d);
    nodalBasis1d->getMassMatrix(tempMassMatrix1d);
    Eigen::MatrixXd massMatrix1d(nlocal1d, nlocal1d);
    copyLuceeToEigen(tempMassMatrix1d, massMatrix1d);

    massMatrixInv1d = massMatrix1d.inverse();
  }

  Lucee::UpdaterStatus
  ProductProjectionCalc1DFrom3D::update(double t)
  {
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    // Distribution function
    const Lucee::Field<3, double>& distfIn = this->getInp<Lucee::Field<3, double> >(0);
    // Field to multiply distribution function by before integrating
    const Lucee::Field<3, double>& distgIn = this->getInp<Lucee::Field<3, double> >(1);
    // Output field
    Lucee::Field<1, double>& resultGlobal = this->getOut<Lucee::Field<1, double> >(0);
    
    // clear out contents of output field
    if (!resultLocal)
    {
      // allocate memory for local moment calculation if not already
      // done: we need to ensure space is also allocated for the
      // ghost-cells as otherwise there is a size mis-match in the
      // allReduce call to sync across velocity space
      Lucee::Region<1, int> localRgn = resultGlobal.getRegion();
      Lucee::Region<1, int> localExtRgn = resultGlobal.getExtRegion();
      
      int lowerConf[1];
      int upperConf[1];
      int lg[1];
      int ug[1];
      for (int i=0; i < 1; ++i)
      {
        lowerConf[i] = localRgn.getLower(i);
        upperConf[i] = localRgn.getUpper(i);
        lg[i] = localRgn.getLower(i) - localExtRgn.getLower(i);
        ug[i] = localExtRgn.getUpper(i) - localRgn.getUpper(i);
      }
      Lucee::Region<1, int> rgnConf(lowerConf, upperConf);
      resultLocal = new Lucee::Field<1, double>(rgnConf, resultGlobal.getNumComponents(), lg, ug);
    }

    // clear out contents of output field
    (*resultLocal) = 0.0;

    Lucee::Region<3, int> globalRgn = grid.getGlobalRegion();
    Lucee::Region<3, int> localRgn = grid.getLocalRegion();
    
    Lucee::ConstFieldPtr<double> distfInPtr = distfIn.createConstPtr();
    Lucee::ConstFieldPtr<double> distgInPtr = distgIn.createConstPtr();
    Lucee::FieldPtr<double> resultPtr = resultLocal->createPtr();

    unsigned nlocal3d = nodalBasis3d->getNumNodes();
    unsigned nlocal1d = nodalBasis1d->getNumNodes();
    int nVolQuad3d = nodalBasis3d->getNumGaussNodes();
    int nVolQuad1d = nodalBasis1d->getNumGaussNodes();

    int localPositionCells = resultGlobal.getExtRegion().getVolume();

    int idx[3];

    Lucee::RowMajorSequencer<3> seq(localRgn);
    
    Eigen::VectorXd distfVec(nlocal3d);
    Eigen::VectorXd distgVec(nlocal3d);

    Eigen::VectorXd distfAtQuad(nVolQuad3d);
    Eigen::VectorXd distgAtQuad(nVolQuad3d);
    
    Eigen::VectorXd rhsIntegrals(nlocal1d);
    Eigen::VectorXd solutionVec(nlocal1d);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      resultLocal->setPtr(resultPtr, idx[0]);
      distfIn.setPtr(distfInPtr, idx);
      distgIn.setPtr(distgInPtr, idx);

      for (int nodeIndex = 0; nodeIndex < nlocal3d; nodeIndex++)
      {
        distfVec(nodeIndex) = distfInPtr[nodeIndex];
        distgVec(nodeIndex) = distgInPtr[nodeIndex];
      }

      // Compute fields at quadrature points
      distfAtQuad = interpMatrix3d*distfVec;
      distgAtQuad = interpMatrix3d*distgVec;
      
      for (int basisIndex = 0; basisIndex < nlocal1d; basisIndex++)
      {
        double integrationResult = 0.0;
        // Compute 3d integration
        for (int quadIndex = 0; quadIndex < nVolQuad3d; quadIndex++)
        {
          integrationResult += gaussWeights3d[quadIndex]*interpMatrix1d(quadIndex % nVolQuad1d, basisIndex)*
            scaleFactor*distfAtQuad(quadIndex)*distgAtQuad(quadIndex);
        }
        rhsIntegrals(basisIndex) = integrationResult;
      }
      // Calculate solution
      solutionVec = massMatrixInv1d*rhsIntegrals;
      
      // Accumulate contribution to local storage
      for (int nodeIndex = 0; nodeIndex < nlocal1d; nodeIndex++)
        resultPtr[nodeIndex] += solutionVec(nodeIndex);
    }

    // Above loop computes moments on local phase-space domain. We need to
    // sum across velocity space to get total momentLocal on configuration
    // space.
    
    // we need to get moment communicator of field as updater's moment
    // communicator is same as its grid's moment communicator. In this
    // case, grid is phase-space grid, which is not what we want.
    TxCommBase *momComm = resultGlobal.getMomComm();
    unsigned xsize = localPositionCells*nlocal1d; // amount to communicate
    momComm->allreduce(xsize, &resultLocal->first(), &resultGlobal.first(), TX_SUM);

    return Lucee::UpdaterStatus();
  }

  void
  ProductProjectionCalc1DFrom3D::declareTypes()
  {
    // Input: distribution function (f)
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: field to multiply with the distribution function (g)
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Output: 1d field containing projection of f*g
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }

  void
  ProductProjectionCalc1DFrom3D::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
