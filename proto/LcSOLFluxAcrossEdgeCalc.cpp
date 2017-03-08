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

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis");
    else
      throw Lucee::Except("SOLFluxAcrossEdgeCalc::readInput: Must specify element to use using 'basis'");
 
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

    // get number of nodes
    unsigned nlocal5d = nodalBasis5d->getNumNodes();

    // get lower and upper surface quadrature data for 5d element
    int nSurfQuad5d = nodalBasis5d->getNumSurfGaussNodes();
    Lucee::Matrix<double> tempSurfQuad5d(nSurfQuad5d, nlocal5d);
    Lucee::Matrix<double> tempSurfCoords5d(nSurfQuad5d, 5);

    surfLowerWeights5d = std::vector<double>(nSurfQuad5d);
    nodalBasis5d->getSurfLowerGaussQuadData(2, tempSurfQuad5d, tempSurfCoords5d, surfLowerWeights5d);

    surfLowerQuad5d = Eigen::MatrixXd(nSurfQuad5d, nlocal5d);
    copyLuceeToEigen(tempSurfQuad5d, surfLowerQuad5d);
    // Undo grid scale of quadrature weights
    double computationalScale = grid.getDx(0)*grid.getDx(1)*grid.getDx(3)*grid.getDx(4);
    for (int quadIndex = 0; quadIndex < surfLowerWeights5d.size(); quadIndex++)
      surfLowerWeights5d[quadIndex] = surfLowerWeights5d[quadIndex]/computationalScale;
    
    surfUpperWeights5d = std::vector<double>(nSurfQuad5d);
    nodalBasis5d->getSurfUpperGaussQuadData(2, tempSurfQuad5d, tempSurfCoords5d, surfUpperWeights5d);
    
    surfUpperQuad5d = Eigen::MatrixXd(nSurfQuad5d, nlocal5d);
    copyLuceeToEigen(tempSurfQuad5d, surfUpperQuad5d);
    // Undo grid scale of quadrature weights
    for (int quadIndex = 0; quadIndex < surfUpperWeights5d.size(); quadIndex++)
      surfUpperWeights5d[quadIndex] = surfUpperWeights5d[quadIndex]/computationalScale;
  }

  Lucee::UpdaterStatus
  SOLFluxAcrossEdgeCalc::update(double t)
  {
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // Distribution function
    const Lucee::Field<5, double>& distfIn = this->getInp<Lucee::Field<5, double> >(0);
    const Lucee::Field<5, double>& bFieldIn = this->getInp<Lucee::Field<5, double> >(1);
    const Lucee::Field<5, double>& hamilDerivIn = this->getInp<Lucee::Field<5, double> >(2);
    // Output dynvector containing total flux on an edge
    Lucee::DynVector<double>& fluxVecOut = this->getOut<Lucee::DynVector<double> >(0);
    
    Lucee::Region<5, int> globalRgn = grid.getGlobalRegion();
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();
    
    Lucee::ConstFieldPtr<double> distfPtr = distfIn.createConstPtr();
    Lucee::ConstFieldPtr<double> bFieldPtr = bFieldIn.createConstPtr();
    Lucee::ConstFieldPtr<double> hamilDerivPtr = hamilDerivIn.createConstPtr();

    unsigned nlocal5d = nodalBasis5d->getNumNodes();
    int nSurfQuad5d = nodalBasis5d->getNumSurfGaussNodes();

    double cellCentroid[5];
    int idx[5];
    int gstIdx[5];

    double localLowerSurfaceFlux = 0.0;
    double localUpperSurfaceFlux = 0.0;

    Eigen::VectorXd distfVec(nlocal5d);
    Eigen::VectorXd bFieldVec(nlocal5d);
    Eigen::VectorXd hamilDerivVec(nlocal5d);

    Eigen::VectorXd distfAtQuad(nSurfQuad5d);
    Eigen::VectorXd bFieldAtQuad(nSurfQuad5d);
    Eigen::VectorXd hamilDerivAtQuad(nSurfQuad5d);
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
        // Get the coordinates of cell center
        grid.setIndex(idx);
        grid.getCentroid(cellCentroid);

        // Only want to calculate outward flux on this surface
        if (cellCentroid[3] < 0.0)
        {
          distfIn.setPtr(distfPtr, idx);
          bFieldIn.setPtr(bFieldPtr, idx);
          hamilDerivIn.setPtr(hamilDerivPtr, idx);

          for (int i = 0; i < nlocal5d; i++)
          {
            distfVec(i) = distfPtr[i];
            bFieldVec(i) = bFieldPtr[i];
            hamilDerivVec(i) = hamilDerivPtr[i];
          }

          // Compute three fields at quadrature points
          distfAtQuad = surfLowerQuad5d*distfVec;
          bFieldAtQuad = surfLowerQuad5d*bFieldVec;
          hamilDerivAtQuad = surfLowerQuad5d*hamilDerivVec;

          for (int quadIndex = 0; quadIndex < nSurfQuad5d; quadIndex++)
            localLowerSurfaceFlux += grid.getSurfArea(2)*surfLowerWeights5d[quadIndex]*distfAtQuad(quadIndex)*
              bFieldAtQuad(quadIndex)*hamilDerivAtQuad(quadIndex);
        }
        else if (integrateGhosts == true)
        {
          hamilDerivIn.setPtr(hamilDerivPtr, idx);
          idx[2] = localRgn.getLower(2)-1;
          // Need to get ghost cell contribution for inward flux on this surface
          distfIn.setPtr(distfPtr, idx);
          bFieldIn.setPtr(bFieldPtr, idx);

          for (int i = 0; i < nlocal5d; i++)
          {
            distfVec(i) = distfPtr[i];
            bFieldVec(i) = bFieldPtr[i];
            hamilDerivVec(i) = hamilDerivPtr[i];
          }

          // Compute three fields at quadrature points
          distfAtQuad = surfUpperQuad5d*distfVec;
          bFieldAtQuad = surfUpperQuad5d*bFieldVec;
          hamilDerivAtQuad = surfLowerQuad5d*hamilDerivVec;

          for (int quadIndex = 0; quadIndex < nSurfQuad5d; quadIndex++)
            localLowerSurfaceFlux += grid.getSurfArea(2)*surfUpperWeights5d[quadIndex]*distfAtQuad(quadIndex)*
              bFieldAtQuad(quadIndex)*hamilDerivAtQuad(quadIndex);
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
        // Get the coordinates of cell center
        grid.setIndex(idx);
        grid.getCentroid(cellCentroid);

        // Only want to calculate outward flux on this surface
        if (cellCentroid[3] > 0.0)
        {
          distfIn.setPtr(distfPtr, idx);
          bFieldIn.setPtr(bFieldPtr, idx);
          hamilDerivIn.setPtr(hamilDerivPtr, idx);

          for (int i = 0; i < nlocal5d; i++)
          {
            distfVec(i) = distfPtr[i];
            bFieldVec(i) = bFieldPtr[i];
            hamilDerivVec(i) = hamilDerivPtr[i];
          }

          // Compute three fields at quadrature points
          distfAtQuad = surfUpperQuad5d*distfVec;
          bFieldAtQuad = surfUpperQuad5d*bFieldVec;
          hamilDerivAtQuad = surfUpperQuad5d*hamilDerivVec;

          for (int quadIndex = 0; quadIndex < nSurfQuad5d; quadIndex++)
            localUpperSurfaceFlux += grid.getSurfArea(2)*surfUpperWeights5d[quadIndex]*distfAtQuad(quadIndex)*
              bFieldAtQuad(quadIndex)*hamilDerivAtQuad(quadIndex);
        }
        else if (integrateGhosts == true)
        {
          hamilDerivIn.setPtr(hamilDerivPtr, idx);
          idx[2] = localRgn.getUpper(2);
          // Need to get ghost cell contribution for inward flux on this surface
          distfIn.setPtr(distfPtr, idx);
          bFieldIn.setPtr(bFieldPtr, idx);

          for (int i = 0; i < nlocal5d; i++)
          {
            distfVec(i) = distfPtr[i];
            bFieldVec(i) = bFieldPtr[i];
            hamilDerivVec(i) = hamilDerivPtr[i];
          }

          // Compute three fields at quadrature points
          distfAtQuad = surfLowerQuad5d*distfVec;
          bFieldAtQuad = surfLowerQuad5d*bFieldVec;
          hamilDerivAtQuad = surfUpperQuad5d*hamilDerivVec;

          for (int quadIndex = 0; quadIndex < nSurfQuad5d; quadIndex++)
            localUpperSurfaceFlux += grid.getSurfArea(2)*surfLowerWeights5d[quadIndex]*distfAtQuad(quadIndex)*
              bFieldAtQuad(quadIndex)*hamilDerivAtQuad(quadIndex);
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
    // Input: 5d magnetic field
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Input: Numerical parallel velocity derivative of hamiltonian
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
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
