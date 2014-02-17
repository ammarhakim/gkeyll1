/**
 * @file	LcElectromagneticZeroVelocitySource.cpp
 *
 * @brief	Takes an input A(z) and computes a CG representation
 * for a distribution function with zero mean velocity
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcElectromagneticZeroVelocitySource.h>
#include <LcMathPhysConstants.h>
#include <math.h> 

namespace Lucee
{
// set id for module system
  const char *ElectromagneticZeroVelocitySource::id = "ElectromagneticZeroVelocitySource";

  ElectromagneticZeroVelocitySource::ElectromagneticZeroVelocitySource()
    : UpdaterIfc()
  {
  }

  void 
  ElectromagneticZeroVelocitySource::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis");
    else
      throw Lucee::Except("ElectromagneticZeroVelocitySource::readInput: Must specify element to use using 'basis'");

    if (tbl.hasNumber("speciesMass"))
      speciesMass = tbl.getNumber("speciesMass");
    else
      throw Lucee::Except("ElectromagneticZeroVelocitySource::readInput: Must specify speciesMass");

    if (tbl.hasNumber("speciesTemp"))
      speciesTemp = tbl.getNumber("speciesTemp");
    else
      throw Lucee::Except("ElectromagneticZeroVelocitySource::readInput: Must specify speciesTemp");

    if (tbl.hasNumber("speciesCharge"))
      speciesCharge = tbl.getNumber("speciesCharge");
    else
      throw Lucee::Except("ElectromagneticZeroVelocitySource::readInput: Must specify speciesCharge");
  }

  void 
  ElectromagneticZeroVelocitySource::initialize()
  {
    Lucee::UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus 
  ElectromagneticZeroVelocitySource::update(double t)
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
        // Set nodal basis index
        nodalBasis->setIndex(idx);
        // Get nodal coordinates into Lucee matrix
        nodalBasis->getNodalCoordinates(nodeCoordsLucee);
        // Copy into Eigen matrix
        copyLuceeToEigen(nodeCoordsLucee, nodeCoords);

        // Get list of nodes owned exclusively
        std::vector<int> ndIds;
        nodalBasis->getExclusiveNodeIndices(ndIds);

        // Copy into output pointer
        for (int componentIndex = 0; componentIndex < ndIds.size(); componentIndex++)
        {
          double pVal = nodeCoords(ndIds[componentIndex], 1);
          double aVal = aPtr[componentIndex];
          double expArg = -(pVal - speciesCharge*aVal)*(pVal - speciesCharge*aVal)/
            (2*speciesMass*speciesTemp*ELEMENTARY_CHARGE);

          double zVal = nodeCoords(ndIds[componentIndex], 0);

          if (fabs(zVal) < 12.5)
            fPtr[componentIndex] = cos(PI*zVal/25.0)*std::exp(expArg);
          else
            fPtr[componentIndex ] = 0.0;
        }
      }
    }

    return Lucee::UpdaterStatus();
  }

  void
  ElectromagneticZeroVelocitySource::declareTypes()
  {
    // takes input A(x) on 2-D grid (CG FIELD)
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    // outputs f(x,p) (CG FIELD)
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }

  void
  ElectromagneticZeroVelocitySource::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}
