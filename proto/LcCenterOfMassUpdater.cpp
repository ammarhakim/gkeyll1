/**
 * @file	LcCenterOfMassUpdater.cpp
 *
 * @brief	Updater to integrate nodal DG/CG field over complete domain
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGlobals.h>
#include <LcCenterOfMassUpdater.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  template <> const char *CenterOfMassUpdater<1>::id = "CenterOfMass1D";
  template <> const char *CenterOfMassUpdater<2>::id = "CenterOfMass2D";
  template <> const char *CenterOfMassUpdater<3>::id = "CenterOfMass3D";

  template <unsigned NDIM>
  CenterOfMassUpdater<NDIM>::CenterOfMassUpdater()
    : Lucee::UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  CenterOfMassUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

// get hold of element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("CenterOfMassUpdater::readInput: Must specify element to use using 'basis'");

    sharedNodes = false;
// check if nodes are shared
    if (tbl.hasBool("shareCommonNodes"))
      sharedNodes = tbl.getBool("shareCommonNodes");

// get background value to subtract off
    backValue = 0.0;
    if (tbl.hasNumber("background"))
      backValue = tbl.getNumber("background");
  }

  template <unsigned NDIM>
  void
  CenterOfMassUpdater<NDIM>::initialize()
  {
    Lucee::UpdaterIfc::initialize();

// get weights for quadrature
    weights.resize(nodalBasis->getNumNodes());
    nodalBasis->getWeights(weights);
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  CenterOfMassUpdater<NDIM>::update(double t)
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
    Lucee::Matrix<double> nodeCoords(nodalBasis->getNumNodes(), 3);

    double localMom[NDIM];
    for (unsigned d=0; d<NDIM; ++d) localMom[d] = 0.0;
    double localAvg = 0;
// loop, performing integration
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      nodalBasis->setIndex(idx);
      nodalBasis->getNodalCoordinates(nodeCoords);

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

// contribution to average
      for (unsigned k=0; k<nlocal; ++k)
        localAvg += weights[k]*(localVals[k]-backValue);
// contribution to moment
      for (unsigned d=0; d<NDIM; ++d)
      {
        for (unsigned k=0; k<nlocal; ++k)
          localMom[d] += weights[k]*nodeCoords(k,d)*(localVals[k]-backValue);
      }
    }

    double volAvg = localAvg;
    double volMom[NDIM];
    for (unsigned d=0; d<NDIM; ++d) volMom[d] = localMom[d];

// get hold of comm pointer to do all parallel messaging
    TxCommBase *comm = this->getComm();
    comm->allreduce(1, &localAvg, &volAvg, TX_SUM);
    comm->allreduce(NDIM, &localMom[0], &volMom[0], TX_SUM);

    std::vector<double> data(NDIM);
    for (unsigned d=0; d<NDIM; ++d)
      data[d] = volMom[d]/volAvg;
    fldInt.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  CenterOfMassUpdater<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }

// instantiations
  template class Lucee::CenterOfMassUpdater<1>;
  template class Lucee::CenterOfMassUpdater<2>;
  template class Lucee::CenterOfMassUpdater<3>;
}
