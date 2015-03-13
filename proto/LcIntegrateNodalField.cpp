/**
 * @file	LcIntegrateNodalField.cpp
 *
 * @brief	Updater to integrate nodal DG/CG field over complete domain
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGlobals.h>
#include <LcIntegrateNodalField.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  template <> const char *IntegrateNodalField<1>::id = "IntegrateNodalField1D";
  template <> const char *IntegrateNodalField<2>::id = "IntegrateNodalField2D";
  template <> const char *IntegrateNodalField<3>::id = "IntegrateNodalField3D";

  template <unsigned NDIM>
  IntegrateNodalField<NDIM>::IntegrateNodalField()
    : Lucee::UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  IntegrateNodalField<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

// get hold of element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("IntegrateNodalField::readInput: Must specify element to use using 'basis'");

    sharedNodes = false;
// check if nodes are shared
    if (tbl.hasBool("shareCommonNodes"))
      sharedNodes = tbl.getBool("shareCommonNodes");
  }

  template <unsigned NDIM>
  void
  IntegrateNodalField<NDIM>::initialize()
  {
    Lucee::UpdaterIfc::initialize();

// get weights for quadrature
    weights.resize(nodalBasis->getNumNodes());
    nodalBasis->getWeights(weights);
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  IntegrateNodalField<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    const Lucee::Field<NDIM, double>& fld = this->getInp<Lucee::Field<NDIM, double> >(0);
    Lucee::DynVector<double>& fldInt = this->getOut<Lucee::DynVector<double> >(0);

    Lucee::ConstFieldPtr<double> fldPtr = fld.createConstPtr();

    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);

    unsigned nlocal = nodalBasis->getNumNodes();
    std::vector<double> localVals(nlocal);
    int idx[NDIM];

    double localInt = 0.0;
// loop, performing integration
    while (seq.step())
    {
      seq.fillWithIndex(idx);
// set index into element basis
      nodalBasis->setIndex(idx);

// extract values at nodes in cell
      if (sharedNodes)
        nodalBasis->extractFromField(fld, localVals);
      else
      {
        fld.setPtr(fldPtr, idx);
// if nodes are not shared simply copy over data
        for (unsigned k=0; k<nlocal; ++k)
          localVals[k] = fldPtr[k];
      }

// perform quadrature
      for (unsigned k=0; k<nlocal; ++k)
        localInt += weights[k]*localVals[k];
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
  IntegrateNodalField<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }

// instantiations
  template class Lucee::IntegrateNodalField<1>;
  template class Lucee::IntegrateNodalField<2>;
  template class Lucee::IntegrateNodalField<3>;
}
