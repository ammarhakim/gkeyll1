/**
 * @file	LcSOLLenardBernsteinScaleCell3DUpdater.cpp
 *
 * @brief	Accumulates correct amount of an input distribution function (assumed to be a diffusion term)
 * to an existing distribution so that the desired amount of energy in each cell is achieved
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLLenardBernsteinScaleCell3DUpdater.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  const char *SOLLenardBernsteinScaleCell3DUpdater::id = "SOLLenardBernsteinScaleCell3D";

  SOLLenardBernsteinScaleCell3DUpdater::SOLLenardBernsteinScaleCell3DUpdater()
  {
  }

  void
  SOLLenardBernsteinScaleCell3DUpdater::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("SOLLenardBernsteinScaleCell3DUpdater::readInput: Must specify element to use using 'basis3d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis1d"))
      nodalBasis1d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis1d");
    else
      throw Lucee::Except("SOLLenardBernsteinScaleCell3DUpdater::readInput: Must specify element to use using 'basis1d'");
  }

  void
  SOLLenardBernsteinScaleCell3DUpdater::initialize()
  {
    UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus
  SOLLenardBernsteinScaleCell3DUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    // Distribution function after diffusion step
    const Lucee::Field<3, double>& distfDiff = this->getInp<Lucee::Field<3, double> >(0);
    // Energy in each cell after drag step
    const Lucee::Field<1, double>& energyDragIn = this->getInp<Lucee::Field<1, double> >(1);
    // Energy in each cell after diffusion step
    const Lucee::Field<1, double>& energyDiffIn = this->getInp<Lucee::Field<1, double> >(2);
    // Distribution function with drag already applied and diffusion to be added on in update()
    Lucee::Field<3, double>& distfOut = this->getOut<Lucee::Field<3, double> >(0);
    
    Lucee::ConstFieldPtr<double> distfDiffPtr   = distfDiff.createConstPtr();
    Lucee::ConstFieldPtr<double> energyDragInPtr = energyDragIn.createConstPtr();
    Lucee::ConstFieldPtr<double> energyDiffInPtr = energyDiffIn.createConstPtr();
    // Remember not to clear distfOut
    Lucee::FieldPtr<double> distfOutPtr = distfOut.createPtr(); // Output pointer

    unsigned nlocal3d = nodalBasis3d->getNumNodes();
    unsigned nlocal1d = nodalBasis1d->getNumNodes();
    double dt = t-this->getCurrTime();

    int idx[3];

    Lucee::Region<3, int> globalRgn = grid.getGlobalRegion();
    Lucee::Region<3, int> localRgn = grid.getLocalRegion();

    Lucee::RowMajorSequencer<3> seq(localRgn);

    while (seq.step())
    {
      seq.fillWithIndex(idx);
      energyDragIn.setPtr(energyDragInPtr, idx[0]);
      energyDiffIn.setPtr(energyDiffInPtr, idx[0]);

      // Compute scale factor to scale diffusion term so that total energy is zero
      double scaleFactor = -energyDragInPtr[0]/energyDiffInPtr[0];
      
      if (std::isinf(scaleFactor))
        continue;

      distfDiff.setPtr(distfDiffPtr, idx);
      distfOut.setPtr(distfOutPtr, idx);
      // Set fOut to give the right energy at this node
      for (int nodeIndex = 0; nodeIndex < nlocal3d; nodeIndex++)
        distfOutPtr[nodeIndex] = distfOutPtr[nodeIndex] + dt*scaleFactor*distfDiffPtr[nodeIndex];
    }

    return Lucee::UpdaterStatus();
  }

  void
  SOLLenardBernsteinScaleCell3DUpdater::declareTypes()
  {
    // Input: Distribution function after diffusion step
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: Energy in each cell after drag step
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // Input: Energy in each cell after drag step
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // Distribution function with drag already applied and diffusion to be added on in update()
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
  }

  void
  SOLLenardBernsteinScaleCell3DUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
