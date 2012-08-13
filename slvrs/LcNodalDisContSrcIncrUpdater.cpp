/**
 * @file	LcNodalDisContSrcIncrUpdater.cpp
 *
 * @brief	Updater to compute increment from source terms for use in DG schemes
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcNodalDisContSrcIncrUpdater.h>
#include <LcStructuredGridBase.h>

namespace Lucee
{
  template <> const char *NodalDisContSrcIncrUpdater<1>::id = "NodalDgSrcIncrement1D";
  template <> const char *NodalDisContSrcIncrUpdater<2>::id = "NodalDgSrcIncrement2D";
  template <> const char *NodalDisContSrcIncrUpdater<3>::id = "NodalDgSrcIncrement3D";

  template <unsigned NDIM>
  NodalDisContSrcIncrUpdater<NDIM>::NodalDisContSrcIncrUpdater()
    : Lucee::UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  NodalDisContSrcIncrUpdater<NDIM>::~NodalDisContSrcIncrUpdater()
  {
  }

  template <unsigned NDIM>
  void
  NodalDisContSrcIncrUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method first
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("NodalDisContSrcIncrUpdater::readInput: Must specify element to use using 'basis'");

// get list of  terms to use on RHS
    if ( !tbl.hasTable("terms") )
      Lucee::Except lce("NodalDisContSrcIncrUpdater::readInput: Must specify terms to to use!");

    Lucee::LuaTable trmTbl = tbl.getTable("terms");
    rhs = trmTbl.template getAllObjects<Lucee::PointSourceIfc>();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  NodalDisContSrcIncrUpdater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    const Lucee::Field<NDIM, double>& qIn = this->getInp<Lucee::Field<NDIM, double> >(0);;
    Lucee::Field<NDIM, double>& qOut = this->getOut<Lucee::Field<NDIM, double> >(0);

    double dt = t-this->getCurrTime();

// number of components in system
    unsigned nlocal = nodalBasis->getNumNodes();
    unsigned ncIn = qIn.getNumComponents()/nlocal;
    unsigned ncOut = qOut.getNumComponents()/nlocal;

    ts.resize(ncOut);
    std::vector<double> src(ncOut);
    Lucee::Matrix<double> nodeCoords(nodalBasis->getNumNodes(), 3);
    Lucee::ConstFieldPtr<double> qInPtr = qIn.createConstPtr();
    Lucee::FieldPtr<double> qOutPtr = qOut.createPtr();

    double xc[3];
    int idx[NDIM];
// loop over each cell in extended region
    Lucee::Region<NDIM, int> localExtRgn = qIn.getExtRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localExtRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      qIn.setPtr(qInPtr, idx);
      qOut.setPtr(qOutPtr, idx);

      nodalBasis->setIndex(idx);
      nodalBasis->getNodalCoordinates(nodeCoords);

      for (unsigned n=0; n<nlocal; ++n)
      {
// copy nodal coordinates over
        for (unsigned d=0; d<3; ++d)
          xc[d] = nodeCoords(n,d);
// compute source
        calcSource(t, xc, &qInPtr[ncIn*n], src);
// add increment into output field
        for (unsigned k=0; k<ncOut; ++k)
          qOutPtr[ncOut*n+k] += dt*src[k];
      }
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  NodalDisContSrcIncrUpdater<NDIM>::declareTypes()
  {
// expect one output field
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void
  NodalDisContSrcIncrUpdater<NDIM>::calcSource(double tm, const double xc[3], const double *inp, std::vector<double>& src)
  {
    unsigned n = src.size();
// zap sources first
    for (unsigned k=0; k<n; ++k)
      src[k] = 0.0;
// compute each source, accumulating it
    for (unsigned i=0; i<rhs.size(); ++i)
    {
      for (unsigned k=0; k<n; ++k)
        ts[k] = 0.0;
      rhs[i]->calcSource(tm, xc, inp, &ts[0]);
      for (unsigned k=0; k<n; ++k)
        src[k] += ts[k];
    }
  }

// instantiations
  template class NodalDisContSrcIncrUpdater<1>;
  template class NodalDisContSrcIncrUpdater<2>;
  template class NodalDisContSrcIncrUpdater<3>;
}
