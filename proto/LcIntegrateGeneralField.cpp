/**
 * @file	LcIntegrateGeneralField.cpp
 *
 * @brief	Updater to integrate a DG field over the entire domain.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGlobals.h>
#include <LcIntegrateGeneralField.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  template <> const char *IntegrateGeneralField<1>::id = "IntegrateGeneralField1D";
  template <> const char *IntegrateGeneralField<2>::id = "IntegrateGeneralField2D";
  template <> const char *IntegrateGeneralField<3>::id = "IntegrateGeneralField3D";

  template <unsigned NDIM>
  IntegrateGeneralField<NDIM>::IntegrateGeneralField()
    : Lucee::UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  IntegrateGeneralField<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

// get hold of element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("IntegrateGeneralField::readInput: Must specify element to use using 'basis'");
  }

  template <unsigned NDIM>
  void
  IntegrateGeneralField<NDIM>::initialize()
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
  IntegrateGeneralField<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    const Lucee::Field<NDIM, double>& fld = this->getInp<Lucee::Field<NDIM, double> >(0);
    Lucee::DynVector<double>& fldInt = this->getOut<Lucee::DynVector<double> >(0);

    Lucee::ConstFieldPtr<double> fldPtr = fld.createConstPtr();

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
      fld.setPtr(fldPtr, idx);

      Eigen::VectorXd fieldVec(nlocal);
      for (int i = 0; i < nlocal; i++)
        fieldVec(i) = fldPtr[i];

      // Interpolate data to quadrature points
      Eigen::VectorXd fieldAtQuad = volQuad.interpMat*fieldVec;

      // perform quadrature
      for (int i = 0; i < nVolQuad; i++)
        localInt += volQuad.weights[i]*fieldAtQuad(i);
    }

    double volInt = localInt;
    // get hold of comm pointer to do all parallel messaging
    TxCommBase *comm = Loki::SingletonHolder<Lucee::Globals>
      ::Instance().comm;
    comm->allreduce(1, &localInt, &volInt, TX_SUM);
    
    std::vector<double> data(1);
    data[0] = volInt;
    fldInt.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  IntegrateGeneralField<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }

  template <unsigned NDIM>
  void
  IntegrateGeneralField<NDIM>::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }

// instantiations
  template class Lucee::IntegrateGeneralField<1>;
  template class Lucee::IntegrateGeneralField<2>;
  template class Lucee::IntegrateGeneralField<3>;
}
