/**
 * @file	LcPredicateUpdater.cpp
 *
 * @brief	Set region based on predicateupdater
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcPredicateUpdater.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <vector>

namespace Lucee
{

  template <> const char *PredicateUpdater<1>::id = "PredicateRegion1D";
  template <> const char *PredicateUpdater<2>::id = "PredicateRegion2D";
  template <> const char *PredicateUpdater<3>::id = "PredicateRegion3D";

  template <unsigned NDIM>
  PredicateUpdater<NDIM>::PredicateUpdater()
    : fnPredRef(-1), fnEvalRef(-1)
  {
  }

  template <unsigned NDIM>
  void
  PredicateUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);
// get predicate function 
    fnPredRef = tbl.getFunctionRef("predicate");
    fnEvalRef = tbl.getFunctionRef("evaluate");
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  PredicateUpdater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    Lucee::Field<NDIM, double>& q = this->getOut<Lucee::Field<NDIM, double> >(0);

    int idx[NDIM];
    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;

    Lucee::FieldPtr<double> ptr = q.createPtr();
    Lucee::Region<NDIM, int> localExtRgn = q.getExtRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localExtRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      q.setPtr(ptr, idx);
// TODO
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  PredicateUpdater<NDIM>::declareTypes()
  {
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }


// instantiations
  template class PredicateUpdater<1>;
  template class PredicateUpdater<2>;
  template class PredicateUpdater<3>;
}
