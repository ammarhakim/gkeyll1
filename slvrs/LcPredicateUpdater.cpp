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

    std::vector<double> res(q.getNumComponents());

    double xc[3];
    int idx[NDIM];

    Lucee::FieldPtr<double> ptr = q.createPtr();
    Lucee::Region<NDIM, int> localExtRgn = q.getExtRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localExtRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
// get centroid coordinate
      grid.setIndex(idx);
      grid.getCentroid(xc);
      if (getValue(t, xc, res))
      {
// inside region, so set field
        q.setPtr(ptr, idx);
        for (unsigned k=0; k<ptr.getNumComponents(); ++k)
          ptr[k] = res[k];
      }
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  bool
  PredicateUpdater<NDIM>::getValue(double tm, const double loc[3], std::vector<double>& res)
  {
    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;
// push function object on stack
    lua_rawgeti(*L, LUA_REGISTRYINDEX, fnPredRef);
// push variables on stack
    for (unsigned i=0; i<3; ++i)
      lua_pushnumber(*L, loc[i]);
    lua_pushnumber(*L, tm);
// call function
    if (lua_pcall(*L, 4, 1, 0) != 0)
    {
      std::string err(lua_tostring(*L, -1));
      lua_pop(*L, 1);
      Lucee::Except lce("PredicateUpdater::getValue: ");
      lce << "Problem evaluating function supplied as 'predicate' ";
      lce << std::endl << "[" << err << "]";
      throw lce;
    }
// fetch results
    if (!lua_isboolean(*L, -1))
      throw Lucee::Except("PredicateUpdater::getValue: Return value not a bool");
    bool isTrue = lua_toboolean(*L, -1);
    lua_pop(*L, 1);

    if (isTrue)
    {
// only evaluate function if we are inside region
// push function object on stack
      lua_rawgeti(*L, LUA_REGISTRYINDEX, fnEvalRef);
// push variables on stack
      for (unsigned i=0; i<3; ++i)
        lua_pushnumber(*L, loc[i]);
      lua_pushnumber(*L, tm);
// call function
      if (lua_pcall(*L, 4, res.size(), 0) != 0)
      {
        std::string err(lua_tostring(*L, -1));
        lua_pop(*L, 1);
        Lucee::Except lce("PredicateUpdater::getValue: ");
        lce << "Problem evaluating function supplied as 'evaluate' ";
        lce << std::endl << "[" << err << "]";
        throw lce;
      }
// fetch results
      for (int i=-res.size(); i<0; ++i)
      {
        if (!lua_isnumber(*L, i))
          throw Lucee::Except("FunctionSource::getSource: Return value not a number");
        res[res.size()+i] = lua_tonumber(*L, i);
      }
      lua_pop(*L, 1);
    }

    return isTrue;
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
