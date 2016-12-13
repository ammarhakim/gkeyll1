/**
 * @file	LcSOLTotalIntegralCalc.cpp
 *
 * @brief	Diagnostic to calculate a integral of a function times a distribution function
 * integrated over the entire phase space
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLTotalIntegralCalc.h>

namespace Lucee
{
  template <> const char *SOLTotalIntegralCalc<1>::id = "SOLTotalIntegralCalc1D";
  template <> const char *SOLTotalIntegralCalc<2>::id = "SOLTotalIntegralCalc2D";
  template <> const char *SOLTotalIntegralCalc<3>::id = "SOLTotalIntegralCalc3D";

  template <unsigned NDIM>
  SOLTotalIntegralCalc<NDIM>::SOLTotalIntegralCalc()
  {
  }

  template <unsigned NDIM>
  void
  SOLTotalIntegralCalc<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("SOLTotalIntegralCalc::readInput: Must specify element to use using 'basis'");
 
    if (tbl.hasNumber("scaleFactor"))
      scaleFactor = tbl.getNumber("scaleFactor");
    else scaleFactor = 1.0;
  }

  template <unsigned NDIM>
  void
  SOLTotalIntegralCalc<NDIM>::initialize()
  {
    UpdaterIfc::initialize();

    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    // get number of nodes
    unsigned nlocal = nodalBasis->getNumNodes();

    // get volume quadrature data
    int nVolQuad = nodalBasis->getNumGaussNodes();
    volWeights = std::vector<double>(nVolQuad);
    Lucee::Matrix<double> tempVolQuad(nVolQuad, nlocal);
    Lucee::Matrix<double> tempVolCoords(nVolQuad, NC);

    nodalBasis->getGaussQuadData(tempVolQuad, tempVolCoords, volWeights);

    volQuad = Eigen::MatrixXd(nVolQuad, nlocal);
    copyLuceeToEigen(tempVolQuad, volQuad);
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  SOLTotalIntegralCalc<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    // Distribution function
    const Lucee::Field<NDIM, double>& distfIn = this->getInp<Lucee::Field<NDIM, double> >(0);
    const Lucee::Field<NDIM, double>& weightIn = this->getInp<Lucee::Field<NDIM, double> >(1);
    // Output dynvector
    Lucee::DynVector<double>& energyVecOut = this->getOut<Lucee::DynVector<double> >(0);
    
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    
    Lucee::ConstFieldPtr<double> distfPtr = distfIn.createConstPtr();
    Lucee::ConstFieldPtr<double> weightPtr = weightIn.createConstPtr();

    unsigned nlocal = nodalBasis->getNumNodes();
    int nVolQuad = nodalBasis->getNumGaussNodes();

    int idx[NDIM];

    double localTotalEnergy = 0.0;

    Eigen::VectorXd distfVec(nlocal);
    Eigen::VectorXd weightVec(nlocal);

    Eigen::VectorXd distfAtQuad(nVolQuad);
    Eigen::VectorXd weightAtQuad(nVolQuad);

    // Create a sequencer to loop over phase space
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      grid.setIndex(idx);

      distfIn.setPtr(distfPtr, idx);
      weightIn.setPtr(weightPtr, idx);

      for (int i = 0; i < nlocal; i++)
      {
        distfVec(i) = distfPtr[i];
        weightVec(i) = weightPtr[i];
      }

      // Compute three fields at quadrature points
      distfAtQuad = volQuad*distfVec;
      weightAtQuad = volQuad*weightVec;

      for (int quadIndex = 0; quadIndex < nVolQuad; quadIndex++)
        localTotalEnergy += volWeights[quadIndex]*distfAtQuad(quadIndex)*
          weightAtQuad(quadIndex);
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

  template <unsigned NDIM>
  void
  SOLTotalIntegralCalc<NDIM>::declareTypes()
  {
    // Input: distribution function
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    // Input: weight function
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    // Ouput dynvector
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }

  template <unsigned NDIM>
  void
  SOLTotalIntegralCalc<NDIM>::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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

  // instantiations
  template class SOLTotalIntegralCalc<1>;
  template class SOLTotalIntegralCalc<2>;
  template class SOLTotalIntegralCalc<3>;
}
