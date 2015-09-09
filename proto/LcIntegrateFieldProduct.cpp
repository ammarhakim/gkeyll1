/**
 * @file	LcIntegrateFieldProduct.cpp
 *
 * @brief	Updater to integrate the product of two DG fields on the same basis set over the entire domain.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGlobals.h>
#include <LcIntegrateFieldProduct.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  template <> const char *IntegrateFieldProduct<1>::id = "IntegrateFieldProduct1D";
  template <> const char *IntegrateFieldProduct<2>::id = "IntegrateFieldProduct2D";
  template <> const char *IntegrateFieldProduct<3>::id = "IntegrateFieldProduct3D";

  template <unsigned NDIM>
  IntegrateFieldProduct<NDIM>::IntegrateFieldProduct()
    : Lucee::UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  IntegrateFieldProduct<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    // get hold of element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("IntegrateFieldProduct::readInput: Must specify element to use using 'basis'");
  }

  template <unsigned NDIM>
  void
  IntegrateFieldProduct<NDIM>::initialize()
  {
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    // local region to update
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    
    int nlocal = nodalBasis->getNumNodes();

    // set index to first location in grid (this is okay as in this
    // updater we are assuming grid is uniform)
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    seq.step();
    int idx[NDIM];
    seq.fillWithIndex(idx);
    nodalBasis->setIndex(idx);

    // get data needed for Gaussian quadrature
    int nVolQuad = nodalBasis->getNumGaussNodes();
    std::vector<double> volWeights(nVolQuad);
    Lucee::Matrix<double> tempVolQuad(nVolQuad, nlocal);
    Lucee::Matrix<double> tempVolCoords(nVolQuad, NC);
    volQuad.reset(nVolQuad, nlocal);

    nodalBasis->getGaussQuadData(tempVolQuad, tempVolCoords, volWeights);
    for (int volIndex = 0; volIndex < nVolQuad; volIndex++)
      volQuad.weights(volIndex) = volWeights[volIndex];
    
    copyLuceeToEigen(tempVolQuad, volQuad.interpMat);
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  IntegrateFieldProduct<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    const Lucee::Field<NDIM, double>& fldOne = this->getInp<Lucee::Field<NDIM, double> >(0);
    const Lucee::Field<NDIM, double>& fldTwo = this->getInp<Lucee::Field<NDIM, double> >(1);
    Lucee::DynVector<double>& fldInt = this->getOut<Lucee::DynVector<double> >(0);

    Lucee::ConstFieldPtr<double> fldOnePtr = fldOne.createConstPtr();
    Lucee::ConstFieldPtr<double> fldTwoPtr = fldTwo.createConstPtr();

    // Only loop over the local region
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);

    unsigned nlocal = nodalBasis->getNumNodes();
    int nVolQuad = nodalBasis->getNumGaussNodes();
    int idx[NDIM];

    double localInt = 0.0;
    // loop, performing integration
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      // set index into element basis
      nodalBasis->setIndex(idx);
      fldOne.setPtr(fldOnePtr, idx);
      fldTwo.setPtr(fldTwoPtr, idx);

      Eigen::VectorXd fieldOneVec(nlocal);
      Eigen::VectorXd fieldTwoVec(nlocal);

      for (int i = 0; i < nlocal; i++)
      {
        fieldOneVec(i) = fldOnePtr[i];
        fieldTwoVec(i) = fldTwoPtr[i];
      }

      // Interpolate data to quadrature points
      Eigen::VectorXd fieldOneAtQuad = volQuad.interpMat*fieldOneVec;
      Eigen::VectorXd fieldTwoAtQuad = volQuad.interpMat*fieldTwoVec;

      // perform quadrature
      for (int i = 0; i < nVolQuad; i++)
        localInt += volQuad.weights[i]*fieldOneAtQuad(i)*fieldTwoAtQuad(i);
    }

    double volInt = localInt;
    // get hold of comm pointer to do all parallel messaging
    TxCommBase *comm = this->getComm();
    comm->allreduce(1, &localInt, &volInt, TX_SUM);
    
    std::vector<double> data(1);
    data[0] = volInt;
    fldInt.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  IntegrateFieldProduct<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }

  template <unsigned NDIM>
  void
  IntegrateFieldProduct<NDIM>::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }

// instantiations
  template class Lucee::IntegrateFieldProduct<1>;
  template class Lucee::IntegrateFieldProduct<2>;
  template class Lucee::IntegrateFieldProduct<3>;
}
