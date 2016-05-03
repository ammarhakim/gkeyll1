/**
 * @file	LcSOLPositivityScaleCellUpdater.cpp
 *
 * @brief	Updater used to adjust energy at each cell by use of a drag term
 * Uses result of LcSOLPositivityDragCellUpdater to add the correct portion to match
 * the desired energy at each node.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLPositivityScaleCellUpdater.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  const char *SOLPositivityScaleCellUpdater::id = "SOLPositivityScaleCellUpdater";

  SOLPositivityScaleCellUpdater::SOLPositivityScaleCellUpdater()
  {
  }

  void
  SOLPositivityScaleCellUpdater::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except("SOLPositivityScaleCellUpdater::readInput: Must specify element to use using 'basis5d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("SOLPositivityScaleCellUpdater::readInput: Must specify element to use using 'basis3d'");
  }

  void
  SOLPositivityScaleCellUpdater::initialize()
  {
    UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus
  SOLPositivityScaleCellUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // Distribution function difference between positivity and drag+positivity steps
    const Lucee::Field<5, double>& distfDelta = this->getInp<Lucee::Field<5, double> >(0);
    // Energy at nodes before positivity (target energy)
    const Lucee::Field<3, double>& energyOrigIn = this->getInp<Lucee::Field<3, double> >(1);
    // Energy at nodes after positivity
    const Lucee::Field<3, double>& energyPosIn = this->getInp<Lucee::Field<3, double> >(2);
    // Energy at nodes after maximum drag
    const Lucee::Field<3, double>& energyDragIn = this->getInp<Lucee::Field<3, double> >(3);
    // Distribution function with positivity and correct drag applied
    Lucee::Field<5, double>& distfOut = this->getOut<Lucee::Field<5, double> >(0);

    Lucee::Region<5, int> globalRgn = grid.getGlobalRegion();
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();
    
    Lucee::ConstFieldPtr<double> distfDeltaPtr   = distfDelta.createConstPtr();
    Lucee::ConstFieldPtr<double> energyOrigInPtr = energyOrigIn.createConstPtr();
    Lucee::ConstFieldPtr<double> energyPosInPtr  = energyPosIn.createConstPtr();
    Lucee::ConstFieldPtr<double> energyDragInPtr = energyDragIn.createConstPtr();

    // Remember not to clear distfOut
    Lucee::FieldPtr<double> distfOutPtr = distfOut.createPtr(); // Output pointer

    unsigned nlocal5d = nodalBasis5d->getNumNodes();
    unsigned nlocal3d = nodalBasis3d->getNumNodes();

    double cellCentroid[5];
    int idx[5];

    Lucee::RowMajorSequencer<5> seq(localRgn);

    while (seq.step())
    {
      seq.fillWithIndex(idx);

      energyOrigIn.setPtr(energyOrigInPtr, idx[0], idx[1], idx[2]);
      energyPosIn.setPtr(energyPosInPtr, idx[0], idx[1], idx[2]);
      energyDragIn.setPtr(energyDragInPtr, idx[0], idx[1], idx[2]);

      if (std::fabs(energyPosInPtr[0] - energyOrigInPtr[0]) > 
        1e-10*energyOrigInPtr[0])
      {
        // Compute scale factor to add delta f at this point
        double scaleFactor = (energyOrigInPtr[0]-energyPosInPtr[0])/
          (energyDragInPtr[0]-energyPosInPtr[0]);
        
        // Make sure scaleFactor is less than one
        if (scaleFactor > 1.0)
        {
          //std::cout << "(" << idx[0] << "," << idx[1] << "," << idx[2] << "," << idx[3] << "," << idx[4] << ","
          //  << ") Drag term is larger than possible = " << scaleFactor << std::endl;
          scaleFactor = 1.0;
        }

        distfDelta.setPtr(distfDeltaPtr, idx);
        distfOut.setPtr(distfOutPtr, idx);
        // Set fOut to give the right energy at this node
        for (int nodeIndex = 0; nodeIndex < nlocal5d; nodeIndex++)
        {
          double startVal = distfOutPtr[nodeIndex];
          distfOutPtr[nodeIndex] = distfOutPtr[nodeIndex] + scaleFactor*distfDeltaPtr[nodeIndex];
          if (distfOutPtr[nodeIndex] < 0.0)
          {
            std::cout << "(" << idx[0] << "," << idx[1] << "," << idx[2] << "," << idx[3] << "," << idx[4] << ","
            << ") distfOutPtr[" << nodeIndex << "] = " << distfOutPtr[nodeIndex] << std::endl <<
            "scaleFactor = " << scaleFactor << ", distf = " << startVal << std::endl <<
            "delta = " << distfDeltaPtr[nodeIndex] << std::endl <<
            "energyOrigInPtr = " << energyOrigInPtr[0] << std::endl <<
            "energyPosInPtr = " << energyPosInPtr[0] << std::endl << 
            "energyDragInPtr = " << energyDragInPtr[0] << std::endl;
            distfOutPtr[nodeIndex] = startVal;
          }
        }
      }
    }

    return Lucee::UpdaterStatus();
  }

  void
  SOLPositivityScaleCellUpdater::declareTypes()
  {
    // Input: difference in distribution functions of positivity and drag steps
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Input: target energy of distribution function at nodes
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: energy of distribition function at nodes after positivity
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: energy of distribition function at nodes after maximum drag
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Output: non-negative distribution function with correct energy
    this->appendOutVarType(typeid(Lucee::Field<5, double>));
  }

  void
  SOLPositivityScaleCellUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
