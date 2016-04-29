/**
 * @file	LcSOLTotalEnergyCalc.cpp
 *
 * @brief	Debug object to investigate number conservation by calculating integrated
 * outward flux on a surface
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLTotalEnergyCalc.h>

namespace Lucee
{
  const char *SOLTotalEnergyCalc::id = "SOLTotalEnergyCalc";

  SOLTotalEnergyCalc::SOLTotalEnergyCalc()
  {
  }

  void
  SOLTotalEnergyCalc::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis");
    else
      throw Lucee::Except("SOLTotalEnergyCalc::readInput: Must specify element to use using 'basis'");
 
    if (tbl.hasNumber("scaleFactor"))
      scaleFactor = tbl.getNumber("scaleFactor");
    else scaleFactor = 1.0;
  }

  void
  SOLTotalEnergyCalc::initialize()
  {
    UpdaterIfc::initialize();

    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // get number of nodes
    unsigned nlocal5d = nodalBasis5d->getNumNodes();

    // get volume quadrature data for 5d element
    int nVolQuad5d = nodalBasis5d->getNumGaussNodes();
    volWeights5d = std::vector<double>(nVolQuad5d);
    Lucee::Matrix<double> tempSurfQuad5d(nVolQuad5d, nlocal5d);
    Lucee::Matrix<double> tempSurfCoords5d(nVolQuad5d, 5);

    nodalBasis5d->getGaussQuadData(tempSurfQuad5d, tempSurfCoords5d, volWeights5d);

    volQuad5d = Eigen::MatrixXd(nVolQuad5d, nlocal5d);
    copyLuceeToEigen(tempSurfQuad5d, volQuad5d);
  }

  Lucee::UpdaterStatus
  SOLTotalEnergyCalc::update(double t)
  {
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // Distribution function
    const Lucee::Field<5, double>& distfIn = this->getInp<Lucee::Field<5, double> >(0);
    const Lucee::Field<5, double>& bFieldIn = this->getInp<Lucee::Field<5, double> >(1);
    const Lucee::Field<5, double>& hamilIn = this->getInp<Lucee::Field<5, double> >(2);
    // Output dynvector containing total energy in domain
    Lucee::DynVector<double>& energyVecOut = this->getOut<Lucee::DynVector<double> >(0);
    
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();
    
    Lucee::ConstFieldPtr<double> distfPtr = distfIn.createConstPtr();
    Lucee::ConstFieldPtr<double> bFieldPtr = bFieldIn.createConstPtr();
    Lucee::ConstFieldPtr<double> hamilPtr = hamilIn.createConstPtr();

    unsigned nlocal5d = nodalBasis5d->getNumNodes();
    int nVolQuad5d = nodalBasis5d->getNumGaussNodes();

    double cellCentroid[5];
    int idx[5];

    double localTotalEnergy = 0.0;

    Eigen::VectorXd distfVec(nlocal5d);
    Eigen::VectorXd bFieldVec(nlocal5d);
    Eigen::VectorXd hamilVec(nlocal5d);

    Eigen::VectorXd distfAtQuad(nVolQuad5d);
    Eigen::VectorXd bFieldAtQuad(nlocal5d);
    Eigen::VectorXd hamilAtQuad(nlocal5d);

    // Create a sequencer to loop over (x,y,v,mu) plane
    Lucee::RowMajorSequencer<5> seq(localRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      // Get the coordinates of cell center
      grid.setIndex(idx);
      grid.getCentroid(cellCentroid);

      // Only want to calculate outward flux on this surface
      distfIn.setPtr(distfPtr, idx);
      bFieldIn.setPtr(bFieldPtr, idx);
      hamilIn.setPtr(hamilPtr, idx);

      for (int i = 0; i < nlocal5d; i++)
      {
        distfVec(i) = distfPtr[i];
        bFieldVec(i) = bFieldPtr[i];
        hamilVec(i) = hamilPtr[i];
      }

      // Compute three fields at quadrature points
      distfAtQuad = volQuad5d*distfVec;
      bFieldAtQuad = volQuad5d*bFieldVec;
      hamilAtQuad = volQuad5d*hamilVec;

      for (int quadIndex = 0; quadIndex < nVolQuad5d; quadIndex++)
        localTotalEnergy += volWeights5d[quadIndex]*distfAtQuad(quadIndex)*
          bFieldAtQuad(quadIndex)*hamilAtQuad(quadIndex);
    }

    double totalEnergy = localTotalEnergy;
    // get hold of comm pointer to do all parallel messaging
    TxCommBase *comm = this->getComm();
    comm->allreduce(1, &localTotalEnergy, &totalEnergy, TX_SUM);

    std::vector<double> data(1);
    data[0] = scaleFactor*totalEnergy;

    energyVecOut.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  void
  SOLTotalEnergyCalc::declareTypes()
  {
    // Input: 5d distribution function
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Input: 5d magnetic field
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Input: Numerical parallel velocity derivative of hamiltonian
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Total energy vs. time
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }

  void
  SOLTotalEnergyCalc::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
