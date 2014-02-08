/**
 * @file	LcElectromagneticZeroVelocitySourceProjection.cpp
 *
 * @brief	Takes an input A(z) and computes a DG representation
 * for a distribution function with zero mean velocity
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcElectromagneticZeroVelocitySourceProjection.h>
#include <LcMathPhysConstants.h>
#include <math.h> 

namespace Lucee
{
// set id for module system
  const char *ElectromagneticZeroVelocitySourceProjection::id = "ElectromagneticZeroVelocitySourceProjection";

  ElectromagneticZeroVelocitySourceProjection::ElectromagneticZeroVelocitySourceProjection()
    : UpdaterIfc()
  {
  }

  void 
  ElectromagneticZeroVelocitySourceProjection::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis");
    else
      throw Lucee::Except("ElectromagneticZeroVelocitySourceProjection::readInput: Must specify element to use using 'basis'");

    if (tbl.hasNumber("speciesMass"))
      speciesMass = tbl.getNumber("speciesMass");
    else
      throw Lucee::Except("ElectromagneticZeroVelocitySourceProjection::readInput: Must specify speciesMass");

    if (tbl.hasNumber("speciesTemp"))
      speciesTemp = tbl.getNumber("speciesTemp");
    else
      throw Lucee::Except("ElectromagneticZeroVelocitySourceProjection::readInput: Must specify speciesTemp");

    if (tbl.hasNumber("speciesCharge"))
      speciesCharge = tbl.getNumber("speciesCharge");
    else
      throw Lucee::Except("ElectromagneticZeroVelocitySourceProjection::readInput: Must specify speciesCharge");
  }

  void 
  ElectromagneticZeroVelocitySourceProjection::initialize()
  {
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();
    // global region to update
    Lucee::Region<2, int> globalRgn = grid.getGlobalRegion();

    Lucee::RowMajorSequencer<2> seq(globalRgn);
    seq.step(); // just to get to first index
    int idx[2];
    seq.fillWithIndex(idx);
    nodalBasis->setIndex(idx);
    
    int nlocal = nodalBasis->getNumNodes();
    int numQuadNodes = nodalBasis->getNumGaussNodes();

    // Get interpolation matrix, gaussian quadrature points, and weights
    Lucee::Matrix<double> interpMatrixLucee(numQuadNodes, nlocal);

    Lucee::Matrix<double> gaussOrdinatesLucee(numQuadNodes, 3);

    gaussWeights = std::vector<double>(numQuadNodes);

    // Allocate Eigen matrices
    interpMatrix = Eigen::MatrixXd(numQuadNodes, nlocal);
    gaussOrdinates = Eigen::MatrixXd(numQuadNodes, 3);

    // Get the interpolation matrix for the volume quadrature points
    nodalBasis->getGaussQuadData(interpMatrixLucee, gaussOrdinatesLucee, gaussWeights);
    copyLuceeToEigen(gaussOrdinatesLucee, gaussOrdinates);

    copyLuceeToEigen(interpMatrixLucee, interpMatrix);

    // Compute and store transpose of interpolation matrix
    interpMatrixTranspose = interpMatrix.transpose();

    // Get mass matrix and then copy to Eigen format
    Lucee::Matrix<double> massMatrixLucee(nlocal, nlocal);

    Eigen::MatrixXd massMatrix(nlocal, nlocal);
    
    nodalBasis->getMassMatrix(massMatrixLucee);
    copyLuceeToEigen(massMatrixLucee, massMatrix);

    // Compute and store inverse of mass matrix
    massMatrixInv = massMatrix.inverse();
  }

  Lucee::UpdaterStatus 
  ElectromagneticZeroVelocitySourceProjection::update(double t)
  {
    const Lucee::StructuredGridBase<2>& grid
      = this->getGrid<Lucee::StructuredGridBase<2> >();

    // Magnetic potential (A parallel)
    const Lucee::Field<2, double>& aIn = this->getInp<Lucee::Field<2, double> >(0);
    // Output maxwellian that should have approximately zero mean velocity
    Lucee::Field<2, double>& fOut = this->getOut<Lucee::Field<2, double> >(0);

    int nlocal = nodalBasis->getNumNodes();

    Lucee::Region<2, int> globalRgn = grid.getGlobalRegion();
    Lucee::Region<2, int> extRgn = fOut.getExtRegion();

    Lucee::ConstFieldPtr<double> aPtr = aIn.createConstPtr();
    Lucee::FieldPtr<double> fPtr = fOut.createPtr();

    fPtr = 0.0;

    int idx[2];
    double cellCentroid[3];

    Lucee::Matrix<double> nodeCoordsLucee(nlocal, 3);
    Eigen::MatrixXd nodeCoords(nlocal, 3);
    
    // Loop over all cells
    for (int ix = extRgn.getLower(0); ix < extRgn.getUpper(0); ix++)
    {
      idx[0] = ix;

      for (int iv = extRgn.getLower(1); iv < extRgn.getUpper(1); iv++)
      {
        // Set inputs
        aIn.setPtr(aPtr, ix, iv);
        // Set outputs
        fOut.setPtr(fPtr, ix, iv);

        // Set grid location
        idx[1] = iv;
        grid.setIndex(idx);
        grid.getCentroid(cellCentroid);
        // Set nodal basis index
        nodalBasis->setIndex(idx);
        // Get nodal coordinates into Lucee matrix
        nodalBasis->getNodalCoordinates(nodeCoordsLucee);
        // Copy into Eigen matrix
        copyLuceeToEigen(nodeCoordsLucee, nodeCoords);

        // Compute A(z,p) at quadrature points in the cell
        Eigen::VectorXd aVec(nlocal);
        for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
          aVec(componentIndex) = aPtr[componentIndex];

        Eigen::VectorXd aAtQuadPoints = interpMatrix*aVec;

        Eigen::VectorXd rhsAtQuadPoints(aAtQuadPoints.rows());

        for (int quadNodeIndex = 0; quadNodeIndex < aAtQuadPoints.rows(); quadNodeIndex++)
        {
          double pVal = cellCentroid[1] + gaussOrdinates(quadNodeIndex,1)*grid.getDx(1)/2.0;
          double zVal = cellCentroid[0] + gaussOrdinates(quadNodeIndex,0)*grid.getDx(0)/2.0;
          double aVal = aAtQuadPoints(quadNodeIndex);
          double expArg = -(pVal - speciesCharge*aVal)*(pVal - speciesCharge*aVal)/
            (2*speciesMass*speciesTemp*ELEMENTARY_CHARGE);

          if (fabs(zVal) < 12.5)
            rhsAtQuadPoints(quadNodeIndex) = cos(PI*zVal/25.0)*gaussWeights[quadNodeIndex]*std::exp(expArg);
          else
            rhsAtQuadPoints(quadNodeIndex) = 0.0;
        }

        Eigen::VectorXd fWeights = massMatrixInv*interpMatrixTranspose*rhsAtQuadPoints;

        // Copy into output pointer
        for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
          fPtr[componentIndex] = fWeights(componentIndex);
      }
    }

    return Lucee::UpdaterStatus();
  }

  void
  ElectromagneticZeroVelocitySourceProjection::declareTypes()
  {
    // takes input A(x) on 2-D grid (DG FIELD)
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    // outputs f(x,p) (DG FIELD)
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }

  void
  ElectromagneticZeroVelocitySourceProjection::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}
