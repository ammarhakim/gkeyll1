/**
 * @file	LcSOLFluxAcrossEdgeCalc3D.cpp
 *
 * @brief	Debug object to investigate number conservation by calculating integrated
 * product of two 3d fields on a surface
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLFluxAcrossEdgeCalc3D.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <limits>
#include <vector>

namespace Lucee
{
  const char *SOLFluxAcrossEdgeCalc3D::id = "SOLFluxAcrossEdgeCalc3D";

  SOLFluxAcrossEdgeCalc3D::SOLFluxAcrossEdgeCalc3D()
  {
  }

  void
  SOLFluxAcrossEdgeCalc3D::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis");
    else
      throw Lucee::Except("SOLFluxAcrossEdgeCalc3D::readInput: Must specify element to use using 'basis'");
 
    if (tbl.hasNumber("scaleFactor"))
      scaleFactor = tbl.getNumber("scaleFactor");
    else scaleFactor = 1.0;
  }

  void
  SOLFluxAcrossEdgeCalc3D::initialize()
  {
    UpdaterIfc::initialize();

    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    // get number of nodes
    unsigned nlocal3d = nodalBasis3d->getNumNodes();

    // get lower and upper surface quadrature data for 3d element
    int nSurfQuad3d = nodalBasis3d->getNumSurfGaussNodes();
    Lucee::Matrix<double> tempSurfQuad3d(nSurfQuad3d, nlocal3d);
    Lucee::Matrix<double> tempSurfCoords3d(nSurfQuad3d, 3);

    surfLowerWeights3d = std::vector<double>(nSurfQuad3d);
    nodalBasis3d->getSurfLowerGaussQuadData(2, tempSurfQuad3d, tempSurfCoords3d, surfLowerWeights3d);

    surfLowerQuad3d = Eigen::MatrixXd(nSurfQuad3d, nlocal3d);
    copyLuceeToEigen(tempSurfQuad3d, surfLowerQuad3d);
    
    surfUpperWeights3d = std::vector<double>(nSurfQuad3d);
    nodalBasis3d->getSurfUpperGaussQuadData(2, tempSurfQuad3d, tempSurfCoords3d, surfUpperWeights3d);
    
    surfUpperQuad3d = Eigen::MatrixXd(nSurfQuad3d, nlocal3d);
    copyLuceeToEigen(tempSurfQuad3d, surfUpperQuad3d);
  }

  Lucee::UpdaterStatus
  SOLFluxAcrossEdgeCalc3D::update(double t)
  {
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    // Two fields to multiply and integrate on each end of the box in direction 2
    const Lucee::Field<3, double>& fieldOneIn = this->getInp<Lucee::Field<3, double> >(0);
    const Lucee::Field<3, double>& fieldTwoIn = this->getInp<Lucee::Field<3, double> >(1);
    // Output dynvector containing integral result on a each edge
    Lucee::DynVector<double>& integralVecOut = this->getOut<Lucee::DynVector<double> >(0);
    
    Lucee::Region<3, int> globalRgn = grid.getGlobalRegion();
    Lucee::Region<3, int> localRgn = grid.getLocalRegion();
    
    Lucee::ConstFieldPtr<double> fieldOnePtr = fieldOneIn.createConstPtr();
    Lucee::ConstFieldPtr<double> fieldTwoPtr = fieldTwoIn.createConstPtr();

    unsigned nlocal3d = nodalBasis3d->getNumNodes();
    int nSurfQuad3d = nodalBasis3d->getNumSurfGaussNodes();

    int idx[3];
    int gstIdx[3];

    double localLowerSurfaceIntegral = 0.0;
    double localUpperSurfaceIntegral = 0.0;

    Eigen::VectorXd fieldOneVec(nlocal3d);
    Eigen::VectorXd fieldTwoVec(nlocal3d);

    Eigen::VectorXd fieldOneAtQuad(nSurfQuad3d);
    Eigen::VectorXd fieldTwoAtQuad(nSurfQuad3d);
    // Check to see if we should integrate on the lower z plane
    if (localRgn.getLower(2) == globalRgn.getLower(2))
    {
      // Create a sequencer to loop over (x,y) plane
      Lucee::RowMajorSequencer<3> seqLowerDim(localRgn.deflate(2));
      while (seqLowerDim.step())
      {
        seqLowerDim.fillWithIndex(idx);
        // Set the deflated index to the proper value
        idx[2] = localRgn.getLower(2);
        // Get the coordinates of cell center
        grid.setIndex(idx);

        fieldOneIn.setPtr(fieldOnePtr, idx);
        fieldTwoIn.setPtr(fieldTwoPtr, idx);

        for (int i = 0; i < nlocal3d; i++)
        {
          fieldOneVec(i) = fieldOnePtr[i];
          fieldTwoVec(i) = fieldTwoPtr[i];
        }

        // Compute three fields at quadrature points
        fieldOneAtQuad = surfLowerQuad3d*fieldOneVec;
        fieldTwoAtQuad = surfLowerQuad3d*fieldTwoVec;

        for (int quadIndex = 0; quadIndex < nSurfQuad3d; quadIndex++)
          localLowerSurfaceIntegral += surfLowerWeights3d[quadIndex]*fieldOneAtQuad(quadIndex)*
            fieldTwoAtQuad(quadIndex);
      }
    }

    // Check to see if we should integrate on the upper z plane
    if (localRgn.getUpper(2) == globalRgn.getUpper(2))
    {
      // Create a sequencer to loop over (x,y,v,mu) plane
      Lucee::RowMajorSequencer<3> seqLowerDim(localRgn.deflate(2));
      while (seqLowerDim.step())
      {
        seqLowerDim.fillWithIndex(idx);
        // Set the deflated index to the proper value
        idx[2] = localRgn.getUpper(2)-1;
        // Get the coordinates of cell center
        grid.setIndex(idx);

        fieldOneIn.setPtr(fieldOnePtr, idx);
        fieldTwoIn.setPtr(fieldTwoPtr, idx);

        for (int i = 0; i < nlocal3d; i++)
        {
          fieldOneVec(i) = fieldOnePtr[i];
          fieldTwoVec(i) = fieldTwoPtr[i];
        }

        // Compute three fields at quadrature points
        fieldOneAtQuad = surfUpperQuad3d*fieldOneVec;
        fieldTwoAtQuad = surfUpperQuad3d*fieldTwoVec;
        
        for (int quadIndex = 0; quadIndex < nSurfQuad3d; quadIndex++)
          localUpperSurfaceIntegral += surfUpperWeights3d[quadIndex]*fieldOneAtQuad(quadIndex)*
            fieldTwoAtQuad(quadIndex);
      }
    }

    double totalLowerSurfaceIntegral = localLowerSurfaceIntegral;
    double totalUpperSurfaceIntegral = localUpperSurfaceIntegral;
    // get hold of comm pointer to do all parallel messaging
    TxCommBase *comm = this->getComm();
    comm->allreduce(1, &localLowerSurfaceIntegral, &totalLowerSurfaceIntegral, TX_SUM);
    comm->allreduce(1, &localUpperSurfaceIntegral, &totalUpperSurfaceIntegral, TX_SUM);

    std::vector<double> data(2);
    data[0] = scaleFactor*totalLowerSurfaceIntegral;
    data[1] = scaleFactor*totalUpperSurfaceIntegral;

    integralVecOut.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  void
  SOLFluxAcrossEdgeCalc3D::declareTypes()
  {
    // Two fields to multiply and integrate on each end of the box in direction 2
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Output dynvector containing integral result on a each edge
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }

  void
  SOLFluxAcrossEdgeCalc3D::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
