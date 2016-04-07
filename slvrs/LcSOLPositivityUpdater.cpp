/**
 * @file	LcPositivityUpdater.cpp
 *
 * @brief	Updater to enforce positivity preservation for 5d SOL simulations.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAlignedRectCoordSys.h>
#include <LcField.h>
#include <LcLinAlgebra.h>
#include <LcMathLib.h>
#include <LcSOLPositivityUpdater.h>
#include <LcStructuredGridBase.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
// set id for module system
  const char *SOLPositivityUpdater::id = "SOLPositivityUpdater";

  SOLPositivityUpdater::SOLPositivityUpdater()
    : UpdaterIfc()
  {
  }

  void 
  SOLPositivityUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    // get hold of 5D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except(
        "SOLPositivityUpdater::readInput: Must specify 5D element to use using 'basis5d'");

    // get hold of 3D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except(
        "SOLPositivityUpdater::readInput: Must specify 3D element to use using 'basis3d'");
  }

  void 
  SOLPositivityUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // local region to update
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();

    // get number of nodes in 3D and 5D
    unsigned nlocal3d = nodalBasis3d->getNumNodes();
    unsigned nlocal5d = nodalBasis5d->getNumNodes();

    // get volume interpolation matrices for 3d element
    int nVolQuad3d = nodalBasis3d->getNumGaussNodes();
    std::vector<double> volWeights3d(nVolQuad3d);
    Lucee::Matrix<double> tempVolQuad3d(nVolQuad3d, nlocal3d);
    Lucee::Matrix<double> tempVolCoords3d(nVolQuad3d, 3);

    nodalBasis3d->getGaussQuadData(tempVolQuad3d, tempVolCoords3d, volWeights3d);

    Eigen::MatrixXd volQuad3d(nVolQuad3d, nlocal3d);
    copyLuceeToEigen(tempVolQuad3d, volQuad3d);

    // get volume interpolation matrices for 5d element
    int nVolQuad5d = nodalBasis5d->getNumGaussNodes();
    std::vector<double> volWeights5d(nVolQuad5d);
    Lucee::Matrix<double> tempVolQuad5d(nVolQuad5d, nlocal5d);
    Lucee::Matrix<double> tempVolCoords5d(nVolQuad5d, 5);

    nodalBasis5d->getGaussQuadData(tempVolQuad5d, tempVolCoords5d, volWeights5d);

    Eigen::MatrixXd volQuad5d(nVolQuad5d, nlocal5d);
    copyLuceeToEigen(tempVolQuad5d, volQuad5d);

    densityMatrix = Eigen::MatrixXd(nlocal5d, nlocal3d);
    // Each row is a the integral of a single 5d basis function times each 3d basis function over the cell
    for (int i = 0; i < nlocal5d; i++)
    {
      for (int j = 0; j < nlocal3d; j++)
      {
        // Compute integral of phi5d_i * phi3d_j
        double integralResult = 0.0;
        for (int gaussIndex = 0; gaussIndex < volWeights5d.size(); gaussIndex++)
        {
          integralResult += volWeights5d[gaussIndex]*volQuad5d(gaussIndex, i)*
            volQuad3d(gaussIndex % nVolQuad3d, j);
        }
        densityMatrix(i, j) = integralResult;
      }
    }
  }

  Lucee::UpdaterStatus 
  SOLPositivityUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<5>& grid
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // Field containing jacobian of integration
    const Lucee::Field<3, double>& jacobianField = this->getInp<Lucee::Field<3, double> >(0);
    // Distribution function to be modified
    Lucee::Field<5, double>& distfOut = this->getOut<Lucee::Field<5, double> >(0);

    int nlocal5d = nodalBasis5d->getNumNodes();
    int nlocal3d = nodalBasis3d->getNumNodes();

    Lucee::ConstFieldPtr<double> jacobianPtr = jacobianField.createConstPtr();
    Lucee::FieldPtr<double> distfPtr = distfOut.createPtr();
    
    // local region to update
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();
    int idx[5];

    Lucee::RowMajorSequencer<5> seq(localRgn);
    
    // Used to sequencer
    Eigen::VectorXd jacobianVector(nlocal3d);
    Eigen::VectorXd distfVector(nlocal5d);

    while (seq.step())
    {
      seq.fillWithIndex(idx);
      grid.setIndex(idx);
      distfOut.setPtr(distfPtr, idx);
      jacobianField.setPtr(jacobianPtr, idx[0], idx[1], idx[2]);

      // Fill out jacobianVector
      for (int i = 0; i < nlocal3d; i++)
        jacobianVector(i) = jacobianPtr[i];

      // Fill out distfVector
      for (int i = 0; i < nlocal5d; i++)
        distfVector(i) = distfPtr[i];

      // Compute density
      double originalNum = distfVector.dot(densityMatrix*jacobianVector);

      // Zero out distfVector entries that are negative
      for (int i = 0; i < nlocal5d; i++)
      {
        if (distfVector(i) < 0.0)
          distfVector(i) = 0.0;
      }

      // Compute modified density. This will be greater than originalNum
      double modifiedNum = distfVector.dot(densityMatrix*jacobianVector);

      // Write modified values to distfOut

      for (int i = 0; i < nlocal5d; i++)
        distfPtr[i] = (originalNum/modifiedNum)*distfVector(i);
    }

    return Lucee::UpdaterStatus();
  }

  void
  SOLPositivityUpdater::declareTypes()
  {
    // takes one input: jacobian field
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // returns one output: modified field (e.g. distribution function)
    this->appendOutVarType(typeid(Lucee::Field<5, double>));
  }

  void
  SOLPositivityUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}
