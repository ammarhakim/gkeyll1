/**
 * @file	LcSOLFluxAcrossEdgeCalc.cpp
 *
 * @brief	Debug object to investigate number conservation by calculating integrated
 * outward flux on a surface
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLFluxAcrossEdgeCalc.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <limits>
#include <vector>

namespace Lucee
{
  const char *SOLFluxAcrossEdgeCalc::id = "SOLFluxAcrossEdgeCalc";

  SOLFluxAcrossEdgeCalc::SOLFluxAcrossEdgeCalc()
  {
  }

  void
  SOLFluxAcrossEdgeCalc::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except("SOLFluxAcrossEdgeCalc::readInput: Must specify element to use using 'basis5d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("SOLFluxAcrossEdgeCalc::readInput: Must specify element to use using 'basis3d'");
 
    if (tbl.hasNumber("scaleFactor"))
      scaleFactor = tbl.getNumber("scaleFactor");
    else scaleFactor = 1.0;

    // Indicates if we will integrate in ghost cells in position space
    integrateGhosts = false;
    if (tbl.hasBool("integrateGhosts"))
      integrateGhosts = tbl.getBool("integrateGhosts");
  }

  void
  SOLFluxAcrossEdgeCalc::initialize()
  {
    UpdaterIfc::initialize();

    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // get number of nodes in 3D and 5D
    unsigned nlocal3d = nodalBasis3d->getNumNodes();
    unsigned nlocal5d = nodalBasis5d->getNumNodes();

    // get lower and upper surface quadrature data for 3d element
    int nSurfQuad3d = nodalBasis3d->getNumSurfGaussNodes();
    std::vector<double> surfLowerWeights3d(nSurfQuad3d);
    Lucee::Matrix<double> tempSurfQuad3d(nSurfQuad3d, nlocal3d);
    Lucee::Matrix<double> tempSurfCoords3d(nSurfQuad3d, 3);

    nodalBasis3d->getSurfLowerGaussQuadData(2, tempSurfQuad3d, tempSurfCoords3d, surfLowerWeights3d);
    Eigen::MatrixXd surfLowerQuad3d(nSurfQuad3d, nlocal3d);
    copyLuceeToEigen(tempSurfQuad3d, surfLowerQuad3d);
    
    std::vector<double> surfUpperWeights3d(nSurfQuad3d);
    nodalBasis3d->getSurfUpperGaussQuadData(2, tempSurfQuad3d, tempSurfCoords3d, surfUpperWeights3d);
    Eigen::MatrixXd surfUpperQuad3d(nSurfQuad3d, nlocal3d);
    copyLuceeToEigen(tempSurfQuad3d, surfUpperQuad3d);

    // get lower and upper surface quadrature data for 5d element
    int nSurfQuad5d = nodalBasis5d->getNumSurfGaussNodes();
    std::vector<double> surfLowerWeights5d(nSurfQuad5d);
    Lucee::Matrix<double> tempSurfQuad5d(nSurfQuad5d, nlocal5d);
    Lucee::Matrix<double> tempSurfCoords5d(nSurfQuad5d, 5);

    nodalBasis5d->getSurfLowerGaussQuadData(2, tempSurfQuad5d, tempSurfCoords5d, surfLowerWeights5d);

    Eigen::MatrixXd surfLowerQuad5d(nSurfQuad5d, nlocal5d);
    copyLuceeToEigen(tempSurfQuad5d, surfLowerQuad5d);
    Eigen::MatrixXd surfLowerCoords5d(nSurfQuad5d, nlocal5d);
    copyLuceeToEigen(tempSurfCoords5d, surfLowerCoords5d);
    
    std::vector<double> surfUpperWeights5d(nSurfQuad5d);
    nodalBasis5d->getSurfUpperGaussQuadData(2, tempSurfQuad5d, tempSurfCoords5d, surfUpperWeights5d);
    Eigen::MatrixXd surfUpperQuad5d(nSurfQuad5d, nlocal5d);
    copyLuceeToEigen(tempSurfQuad5d, surfUpperQuad5d);
    Eigen::MatrixXd surfUpperCoords5d(nSurfQuad5d, nlocal5d);
    copyLuceeToEigen(tempSurfCoords5d, surfUpperCoords5d);

    mom0MatrixLower = Eigen::MatrixXd(nlocal5d, nlocal3d);
    mom0MatrixUpper = Eigen::MatrixXd(nlocal5d, nlocal3d);
    mom1MatrixLower = Eigen::MatrixXd(nlocal5d, nlocal3d);
    mom1MatrixUpper = Eigen::MatrixXd(nlocal5d, nlocal3d);
    // Each row is a the integral of a single 5d basis function times each 3d basis function over the cell
    for (int i = 0; i < nlocal5d; i++)
    {
      for (int j = 0; j < nlocal3d; j++)
      {
        // Compute integral of phi5d_i * phi3d_j on lower and upper surfaces
        double mom0ResultLower = 0.0;
        double mom0ResultUpper = 0.0;
        double mom1ResultLower = 0.0;
        double mom1ResultUpper = 0.0;

        // Loop over each 5d quadrature point on the fixed z surface
        for (int gaussIndex = 0; gaussIndex < surfLowerWeights5d.size(); gaussIndex++)
        {
          mom0ResultLower += surfLowerWeights5d[gaussIndex]*surfLowerQuad5d(gaussIndex, i)*
            surfLowerQuad3d(gaussIndex % nSurfQuad3d, j);
          mom0ResultUpper += surfUpperWeights5d[gaussIndex]*surfUpperQuad5d(gaussIndex, i)*
            surfUpperQuad3d(gaussIndex % nSurfQuad3d, j);

          // Parallel velocity at this quadrature point
          double vCoordLower = surfLowerCoords5d(gaussIndex,3)*0.5*grid.getDx(3);
          double vCoordUpper = surfUpperCoords5d(gaussIndex,3)*0.5*grid.getDx(3);

          mom1ResultLower += surfLowerWeights5d[gaussIndex]*surfLowerQuad5d(gaussIndex, i)*
            surfLowerQuad3d(gaussIndex % nSurfQuad3d, j)*vCoordLower;
          mom1ResultUpper += surfUpperWeights5d[gaussIndex]*surfUpperQuad5d(gaussIndex, i)*
            surfUpperQuad3d(gaussIndex % nSurfQuad3d, j)*vCoordUpper;
        }
        mom0MatrixLower(i, j) = mom0ResultLower;
        mom0MatrixUpper(i, j) = mom0ResultUpper;
        mom1MatrixLower(i, j) = mom1ResultLower;
        mom1MatrixUpper(i, j) = mom1ResultUpper;
      }
    }
  }

  Lucee::UpdaterStatus
  SOLFluxAcrossEdgeCalc::update(double t)
  {
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // Distribution function
    const Lucee::Field<5, double>& distfIn = this->getInp<Lucee::Field<5, double> >(0);
    const Lucee::Field<3, double>& bFieldIn = this->getInp<Lucee::Field<3, double> >(1);
    // Output dynvector containing total flux on an edge
    Lucee::DynVector<double>& fluxVecOut = this->getOut<Lucee::DynVector<double> >(0);
    
    Lucee::Region<5, int> globalRgn = grid.getGlobalRegion();
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();
    
    Lucee::ConstFieldPtr<double> distfPtr = distfIn.createConstPtr();
    Lucee::ConstFieldPtr<double> bFieldPtr = bFieldIn.createConstPtr();

    unsigned nlocal3d = nodalBasis3d->getNumNodes();
    unsigned nlocal5d = nodalBasis5d->getNumNodes();

    double cellCentroid[5];
    int idx[5];
    int sknIdx[5];
    int gstIdx[5];

    double localLowerSurfaceFlux = 0.0;
    double localUpperSurfaceFlux = 0.0;

    // Check to see if we should integrate on the lower z plane
    if (localRgn.getLower(2) == globalRgn.getLower(2))
    {
      // Create a sequencer to loop over (x,y,v,mu) plane
      Lucee::RowMajorSequencer<5> seqLowerDim(localRgn.deflate(2));
      while (seqLowerDim.step())
      {
        seqLowerDim.fillWithIndex(idx);
        // Set the deflated index to the proper value
        idx[2] = localRgn.getLower(2);
        // Set skin and ghost cell indices
        for (int i = 0; i < 5; i++)
        {
          sknIdx[i] = idx[i];
          gstIdx[i] = idx[i];
        }
        gstIdx[2] = idx[2] - 1;
        // Get the coordinates of cell center
        grid.setIndex(idx);
        grid.getCentroid(cellCentroid);
        // Fill out solution and magnetic field at this location
        distfIn.setPtr(distfPtr, idx);
        bFieldIn.setPtr(bFieldPtr, idx[0], idx[1], idx[2]);
        Eigen::VectorXd distfVec(nlocal5d);
        Eigen::VectorXd bFieldVec(nlocal3d);
        for (int i = 0; i < nlocal5d; i++)
          distfVec(i) = distfPtr[i];
        for (int i = 0; i < nlocal3d; i++)
          bFieldVec(i) = bFieldPtr[i];

        // Only want to calculate outward flux on this surface
        if (cellCentroid[3] < 0.0)
          localLowerSurfaceFlux += distfVec.dot((cellCentroid[3]*mom0MatrixLower + mom1MatrixLower)*bFieldVec);
        else if (integrateGhosts == true)
        {
          // Need to get ghost cell contribution for inward flux on this surface
          distfIn.setPtr(distfPtr, gstIdx);
          for (int i = 0; i < nlocal5d; i++)
            distfVec(i) = distfPtr[i];
          bFieldIn.setPtr(bFieldPtr, gstIdx[0], gstIdx[1], gstIdx[2]);
          for (int i = 0; i < nlocal3d; i++)
            bFieldVec(i) = bFieldPtr[i];
          localLowerSurfaceFlux += distfVec.dot((cellCentroid[3]*mom0MatrixUpper + mom1MatrixUpper)*bFieldVec);
        }
      }
    }

    // Check to see if we should integrate on the upper z plane
    if (localRgn.getUpper(2) == globalRgn.getUpper(2))
    {
      // Create a sequencer to loop over (x,y,v,mu) plane
      Lucee::RowMajorSequencer<5> seqLowerDim(localRgn.deflate(2));
      while (seqLowerDim.step())
      {
        seqLowerDim.fillWithIndex(idx);
        // Set the deflated index to the proper value
        idx[2] = localRgn.getUpper(2)-1;
        // Set skin and ghost cell indices
        for (int i = 0; i < 5; i++)
        {
          sknIdx[i] = idx[i];
          gstIdx[i] = idx[i];
        }
        gstIdx[2] = idx[2] + 1;
        // Get the coordinates of cell center
        grid.setIndex(idx);
        grid.getCentroid(cellCentroid);
        // Fill out solution and magnetic field at this location
        distfIn.setPtr(distfPtr, idx);
        bFieldIn.setPtr(bFieldPtr, idx[0], idx[1], idx[2]);
        Eigen::VectorXd distfVec(nlocal5d);
        Eigen::VectorXd bFieldVec(nlocal3d);
        for (int i = 0; i < nlocal5d; i++)
          distfVec(i) = distfPtr[i];
        for (int i = 0; i < nlocal3d; i++)
          bFieldVec(i) = bFieldPtr[i];

        // Only want to calculate outward flux on this surface
        if (cellCentroid[3] > 0.0)
          localUpperSurfaceFlux += distfVec.dot((cellCentroid[3]*mom0MatrixUpper + mom1MatrixUpper)*bFieldVec);
        else if (integrateGhosts == true)
        {
          // Need to get ghost cell contribution for inward flux on this surface
          distfIn.setPtr(distfPtr, gstIdx);
          for (int i = 0; i < nlocal5d; i++)
            distfVec(i) = distfPtr[i];
          bFieldIn.setPtr(bFieldPtr, gstIdx[0], gstIdx[1], gstIdx[2]);
          for (int i = 0; i < nlocal3d; i++)
            bFieldVec(i) = bFieldPtr[i];
          localUpperSurfaceFlux += distfVec.dot((cellCentroid[3]*mom0MatrixLower + mom1MatrixLower)*bFieldVec);
        }
      }
    }

    double totalLowerSurfaceFlux = localLowerSurfaceFlux;
    double totalUpperSurfaceFlux = localUpperSurfaceFlux;
    // get hold of comm pointer to do all parallel messaging
    TxCommBase *comm = this->getComm();
    comm->allreduce(1, &localLowerSurfaceFlux, &totalLowerSurfaceFlux, TX_SUM);
    comm->allreduce(1, &localUpperSurfaceFlux, &totalUpperSurfaceFlux, TX_SUM);

    std::vector<double> data(2);
    data[0] = scaleFactor*totalLowerSurfaceFlux;
    data[1] = scaleFactor*totalUpperSurfaceFlux;

    fluxVecOut.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  void
  SOLFluxAcrossEdgeCalc::declareTypes()
  {
    // Input: 5d distribution function
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Input: 3d magnetic field
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Total moment on lower and upper surfaces
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }

  void
  SOLFluxAcrossEdgeCalc::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
