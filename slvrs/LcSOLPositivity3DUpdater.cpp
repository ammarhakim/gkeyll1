/**
 * @file	LcPositivityUpdater.cpp
 *
 * @brief	Updater to enforce positivity preservation for 3d SOL simulations.
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
#include <LcSOLPositivity3DUpdater.h>
#include <LcStructuredGridBase.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
// set id for module system
  const char *SOLPositivity3DUpdater::id = "SOLPositivity3DUpdater";

  SOLPositivity3DUpdater::SOLPositivity3DUpdater()
    : UpdaterIfc()
  {
  }

  void 
  SOLPositivity3DUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    // get hold of 3D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis");
    else
      throw Lucee::Except(
        "SOLPositivity3DUpdater::readInput: Must specify 3D element to use using 'basis'");
  }

  void 
  SOLPositivity3DUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();
    
    int nlocal = nodalBasis->getNumNodes();
    // Get 3D Mass Matrix
    Lucee::Matrix<double> tempMassMatrix3d(nlocal, nlocal);
    nodalBasis->getMassMatrix(tempMassMatrix3d);
    Eigen::MatrixXd massMatrix3d(nlocal, nlocal);
    copyLuceeToEigen(tempMassMatrix3d, massMatrix3d);
    // Compute mom0Vector
    mom0Vector = massMatrix3d.colwise().sum();
  }

  Lucee::UpdaterStatus 
  SOLPositivity3DUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<3>& grid
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    // Distribution function to be modified
    Lucee::Field<3, double>& distfOut = this->getOut<Lucee::Field<3, double> >(0);

    int nlocal = nodalBasis->getNumNodes();

    Lucee::FieldPtr<double> distfPtr = distfOut.createPtr();
    
    // local region to update
    Lucee::Region<3, int> localRgn = grid.getLocalRegion();
    int idx[3];

    Lucee::RowMajorSequencer<3> seq(localRgn);
    
    // Used in sequencer
    Eigen::VectorXd distfVector(nlocal);

    while (seq.step())
    {
      seq.fillWithIndex(idx);
      grid.setIndex(idx);
      distfOut.setPtr(distfPtr, idx);

      // Fill out distfVector
      for (int i = 0; i < nlocal; i++)
        distfVector(i) = distfPtr[i];

      // Compute density
      double originalNum = mom0Vector.dot(distfVector);
      if (originalNum < 0.0)
      {
        std::cout << "(" << idx[0] << "," << idx[1] << "," << idx[2] << 
          ") entire cell negative (density = " << originalNum << ")" << std::endl;
        //std::cout << distfVector << std::endl;
        //return Lucee::UpdaterStatus(false, 0.0);
        // Write modified values to distfOut
        for (int i = 0; i < nlocal; i++)
          distfPtr[i] = 0.0;
        continue;
      }
      else if (originalNum == 0.0)
      {
        //std::cout << "(positivity) cell is zero. skipping." << std::endl;
        continue;
      }

      // Zero out distfVector entries that are negative
      for (int i = 0; i < nlocal; i++)
      {
        if (distfVector(i) < 0.0)
          distfVector(i) = 0.0;
      }

      // Compute modified density. This will be greater than originalNum
      double modifiedNum = mom0Vector.dot(distfVector);

      if (modifiedNum == 0.0)
      {
        std::cout << "New cell is all zero" << std::endl;
        for (int i = 0; i < nlocal; i++)
          distfPtr[i] = 0.0;
      }
      else
      {
        // Write modified values to distfOut
        for (int i = 0; i < nlocal; i++)
          distfPtr[i] = (originalNum/modifiedNum)*distfVector(i);
      }
    }

    return Lucee::UpdaterStatus();
  }

  void
  SOLPositivity3DUpdater::declareTypes()
  {
    // returns one output: modified field (e.g. distribution function)
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
  }

  void
  SOLPositivity3DUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}
