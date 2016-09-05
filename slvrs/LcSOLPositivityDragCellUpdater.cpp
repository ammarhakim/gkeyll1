/**
 * @file	LcSOLPositivityDragCellUpdater.cpp
 *
 * @brief	Updater used to adjust energy at each cell by use of a drag term
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLPositivityDragCellUpdater.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  const char *SOLPositivityDragCellUpdater::id = "SOLPositivityDragCellUpdater";

  SOLPositivityDragCellUpdater::SOLPositivityDragCellUpdater()
  {
  }

  void
  SOLPositivityDragCellUpdater::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except("SOLPositivityDragCellUpdater::readInput: Must specify element to use using 'basis5d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("SOLPositivityDragCellUpdater::readInput: Must specify element to use using 'basis3d'");

    if (tbl.hasNumber("scaleFactor"))
      scaleFactor = tbl.getNumber("scaleFactor");
    else scaleFactor = 1.0;
  }

  void
  SOLPositivityDragCellUpdater::initialize()
  {
    UpdaterIfc::initialize();

    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

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

    // Get 3D Mass Matrix
    Lucee::Matrix<double> tempMassMatrix3d(nlocal3d, nlocal3d);
    nodalBasis3d->getMassMatrix(tempMassMatrix3d);
    Eigen::MatrixXd massMatrix3d(nlocal3d, nlocal3d);
    copyLuceeToEigen(tempMassMatrix3d, massMatrix3d);
    // Compute mom0Vector
    mom0Vector = massMatrix3d.colwise().sum();

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
  SOLPositivityDragCellUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // Distribution function with incorrect energy
    const Lucee::Field<5, double>& distfIn = this->getInp<Lucee::Field<5, double> >(0);
    // 3d magnetic field (to compute total number in a cell)
    const Lucee::Field<3, double>& bFieldIn = this->getInp<Lucee::Field<3, double> >(1);
    // Energy of distfIn at each node
    const Lucee::Field<3, double>& energyModIn = this->getInp<Lucee::Field<3, double> >(2);
    // Desired energy at each node
    const Lucee::Field<3, double>& energyOrigIn = this->getInp<Lucee::Field<3, double> >(3);
    // Mean velocity field n*u_parallel
    const Lucee::Field<3, double>& nTimesUIn = this->getInp<Lucee::Field<3, double> >(4);
    // Density field
    const Lucee::Field<3, double>& nIn = this->getInp<Lucee::Field<3, double> >(5);
    // Distribution function after drag term
    Lucee::Field<5, double>& distfOut = this->getOut<Lucee::Field<5, double> >(0);

    Lucee::Region<5, int> globalRgn = grid.getGlobalRegion();
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();
    
    Lucee::ConstFieldPtr<double> distfInPtr       = distfIn.createConstPtr();
    // Another pointer to identify upwind neighbor
    Lucee::ConstFieldPtr<double> distfInUpwindPtr = distfIn.createConstPtr();
    Lucee::ConstFieldPtr<double> bFieldInPtr      = bFieldIn.createConstPtr();
    Lucee::ConstFieldPtr<double> energyModInPtr   = energyModIn.createConstPtr();
    Lucee::ConstFieldPtr<double> energyOrigInPtr  = energyOrigIn.createConstPtr();
    Lucee::ConstFieldPtr<double> nTimesUInPtr     = nTimesUIn.createConstPtr();
    Lucee::ConstFieldPtr<double> nInPtr           = nIn.createConstPtr();

    Lucee::FieldPtr<double> distfOutPtr = distfOut.createPtr(); // Output pointer
    distfOut = 0.0; // clear out current contents

    unsigned nlocal5d = nodalBasis5d->getNumNodes();
    unsigned nlocal3d = nodalBasis3d->getNumNodes();

    double cellCentroid[5];
    int idx[5];
    int idxUpwind[5];
    double dt = t-this->getCurrTime();
    // Figure out value of parallel velocity on last edge that flux is computed on
    idx[0] = globalRgn.getUpper(0)-1;
    idx[1] = globalRgn.getUpper(1)-1;
    idx[2] = globalRgn.getUpper(2)-1;
    idx[3] = globalRgn.getUpper(3)-1;
    idx[4] = globalRgn.getUpper(4)-1;
    // Set grid to last parallel velocity cell (upper boundary)
    grid.setIndex(idx);
    grid.getCentroid(cellCentroid);
    double maxPosV = (cellCentroid[3] - 0.5*grid.getDx(3));

    // Used in sequencer
    Eigen::VectorXd bFieldVec(nlocal3d);
    Eigen::VectorXd distfVec(nlocal5d);
    Eigen::VectorXd nVec(nlocal3d);
    Eigen::VectorXd nTimesUVec(nlocal3d);

    Lucee::RowMajorSequencer<5> seq(localRgn);
    
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      // Get the coordinates of cell center
      grid.setIndex(idx);
      grid.getCentroid(cellCentroid);

      energyModIn.setPtr(energyModInPtr, idx[0], idx[1], idx[2]);
      energyOrigIn.setPtr(energyOrigInPtr, idx[0], idx[1], idx[2]);

      // check to see if energy at this cell needs to be corrected
      if (std::fabs(energyModInPtr[0] - energyOrigInPtr[0]) > 
                1e-10*energyOrigInPtr[0])
      {
        bFieldIn.setPtr(bFieldInPtr, idx[0], idx[1], idx[2]);
        nTimesUIn.setPtr(nTimesUInPtr, idx[0], idx[1], idx[2]);
        nIn.setPtr(nInPtr, idx[0], idx[1], idx[2]);

        for (int nodeIndex = 0; nodeIndex < nlocal3d; nodeIndex++)
        {
          bFieldVec(nodeIndex) = bFieldInPtr[nodeIndex];
          nTimesUVec(nodeIndex) = nTimesUInPtr[nodeIndex];
          nVec(nodeIndex) = nInPtr[nodeIndex];
        }

        // Compute average value of u_parallel
        double avgU = (mom0Vector.dot(nTimesUVec))/(mom0Vector.dot(nVec));

        distfIn.setPtr(distfInPtr, idx);
        for (int nodeIndex = 0; nodeIndex < nlocal5d; nodeIndex++)
          distfVec(nodeIndex) = distfInPtr[nodeIndex];

        double fOldAvg = distfVec.dot(densityMatrix*bFieldVec);
        double fLowerAvg = 0.0;
        double fUpperAvg = 0.0;

        if (fOldAvg == 0.0)
        {
          //std::cout << "fOldAvg is zero. This should probably not happen." << std::endl;
          //std::cout << "distfReduced" << std::endl << distfReduced << std::endl;
          continue;
        }
        else if (fOldAvg < 0.0)
        {
          std::cout << "fOldAvg is less than zero. This should probably not happen." << std::endl;
          //std::cout << "distfReduced" << std::endl << distfReduced << std::endl;
        }
        
        // Copy idx to idxUpwind
        for (int dimIndex = 0; dimIndex < 5; dimIndex++)
          idxUpwind[dimIndex] = idx[dimIndex];

        // Figure out which neighboring cell to use to calculate upwindAvg
        if (cellCentroid[3]-avgU > 0.0)
        {
          fLowerAvg = fOldAvg;
          if (idx[3] < globalRgn.getUpper(3)-1)
          {
            idxUpwind[3] = idx[3] + 1;
            distfIn.setPtr(distfInUpwindPtr, idxUpwind);
            for (int nodeIndex = 0; nodeIndex < nlocal5d; nodeIndex++)
              distfVec(nodeIndex) = distfInUpwindPtr[nodeIndex];
            fUpperAvg = distfVec.dot(densityMatrix*bFieldVec);
          }
        }
        else if (cellCentroid[3]-avgU < 0.0)
        {
          fUpperAvg = fOldAvg;
          if (idx[3] > globalRgn.getLower(3))
          {
            idxUpwind[3] = idx[3] - 1;
            distfIn.setPtr(distfInUpwindPtr, idxUpwind);
            for (int nodeIndex = 0; nodeIndex < nlocal5d; nodeIndex++)
              distfVec(nodeIndex) = distfInUpwindPtr[nodeIndex];
            fLowerAvg = distfVec.dot(densityMatrix*bFieldVec);
          }
        }

        // Put in limiters on cell average?
        if (fUpperAvg < 0.0)
          fUpperAvg = 0.0;
        if (fLowerAvg < 0.0)
          fLowerAvg = 0.0;

        double alpha = 1/(maxPosV + std::fabs(avgU));
        double fIncrement = (cellCentroid[3] + 0.5*grid.getDx(3) - avgU)*fUpperAvg -
          (cellCentroid[3] - 0.5*grid.getDx(3) - avgU)*fLowerAvg;
        // Compute new density of cell after maximum drag step
        double fNewAvg = fOldAvg + alpha*fIncrement;

        if (fNewAvg < 0.0)
        {
          // Set fNewAvg to 0 if negative (presumably it happens when answer is approximately zero)
          if (std::fabs(fNewAvg/fOldAvg) < 1e-8)
          {
            fNewAvg = 0.0;
          }
          else
          {
          std::cout << "fNewAvg is negative! = " << fNewAvg << std::endl;
          std::cout << "fOldAvg = " << fOldAvg << std::endl;
          std::cout << "alpha*fIncrement = " << alpha*fIncrement << std::endl;
          std::cout << "fIncrement upper = " << alpha*(cellCentroid[3] + 0.5*grid.getDx(3))*fUpperAvg << std::endl;
          std::cout << "fIncrement lower = " << alpha*(cellCentroid[3] - 0.5*grid.getDx(3))*fLowerAvg << std::endl;
          }
        }

            
        distfOut.setPtr(distfOutPtr, idx);
        // Set fOut by scaling fIn to have the right density
        for (int nodeIndex = 0; nodeIndex < nlocal5d; nodeIndex++)
          distfOutPtr[nodeIndex] = (fNewAvg/fOldAvg)*distfInPtr[nodeIndex];
      }
    }
   
    return Lucee::UpdaterStatus();
  }

  void
  SOLPositivityDragCellUpdater::declareTypes()
  {
    // Input: modified distribution function
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Input: bField3d
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: energy of modified distribution function at cells
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: energy of original distribition function at cells
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: n*u_parallel
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: n
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Output: another modified distribution function
    this->appendOutVarType(typeid(Lucee::Field<5, double>));
  }

  void
  SOLPositivityDragCellUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
