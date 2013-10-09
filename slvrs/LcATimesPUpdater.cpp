/**
 * @file	LcATimesPUpdater.cpp
 *
 * @brief	Projects A(x)^2 onto the same basis functions that A(x) uses
 * Inputs and outputs assume CG fields, not DG fields.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcATimesPUpdater.h>

namespace Lucee
{
// set id for module system
  const char *ATimesPUpdater::id = "ATimesPUpdater";

  ATimesPUpdater::ATimesPUpdater()
    : UpdaterIfc()
  {
  }

  void 
  ATimesPUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis");
    else
      throw Lucee::Except("ATimesPUpdater::readInput: Must specify element to use using 'basis'");
  }

  void 
  ATimesPUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus 
  ATimesPUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<2>& grid
      = this->getGrid<Lucee::StructuredGridBase<2> >();

    // Electron and ion densities
    const Lucee::Field<2, double>& aIn = this->getInp<Lucee::Field<2, double> >(0);
    Lucee::Field<2, double>& aTimesPOut = this->getOut<Lucee::Field<2, double> >(0);

    int nlocal = nodalBasis->getNumNodes();

    Lucee::Region<2, int> globalRgn = grid.getGlobalRegion();

    Lucee::ConstFieldPtr<double> aPtr = aIn.createConstPtr();
    Lucee::FieldPtr<double> aTimesPPtr = aTimesPOut.createPtr();

    aTimesPPtr = 0.0;

    int idx[2];

    Lucee::Matrix<double> nodeCoordsLucee(nlocal, 3);
    Eigen::MatrixXd nodeCoords(nlocal, 3);
    
    // Loop over all cells
    for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
    {
      idx[0] = ix;

      for (int iv = globalRgn.getLower(1); iv < globalRgn.getUpper(1); iv++)
      {
        // Set inputs
        aIn.setPtr(aPtr, ix, iv);
        // Set outputs
        aTimesPOut.setPtr(aTimesPPtr, ix, iv);

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
          aTimesPPtr[componentIndex] = aPtr[componentIndex]*nodeCoords(ndIds[componentIndex], 1);
      }
    }

    return Lucee::UpdaterStatus();
  }

  void
  ATimesPUpdater::declareTypes()
  {
    // takes input A(x) on 2-D grid (CG FIELD)
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    // outputs f(x,p) = A(x)*p (CG FIELD)
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }

  void
  ATimesPUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}
