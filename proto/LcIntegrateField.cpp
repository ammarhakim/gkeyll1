/**
 * @file	LcIntegrateField.cpp
 *
 * @brief	Updater to integrate a field along a specified line.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGlobals.h>
#include <LcIntegrateField.h>
#include <LcLuaState.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <vector>

namespace Lucee
{
  template <> const char *IntegrateField<1>::id = "IntegrateField1D";
  template <> const char *IntegrateField<2>::id = "IntegrateField2D";
  template <> const char *IntegrateField<3>::id = "IntegrateField3D";

  template <unsigned NDIM>
  IntegrateField<NDIM>::IntegrateField()
    : Lucee::UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  IntegrateField<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    hasFunction = false;
// get reference to function
    if (tbl.hasFunction("integrand"))
    {
      hasFunction = true;
      fnRef = tbl.getFunctionRef("integrand");
    }
  }

  template <unsigned NDIM>
  void
  IntegrateField<NDIM>::initialize()
  {
    Lucee::UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  IntegrateField<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    const Lucee::Field<NDIM, double>& fld = this->getInp<Lucee::Field<NDIM, double> >(0);
    Lucee::DynVector<double>& fldInt = this->getOut<Lucee::DynVector<double> >(0);

    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();    
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);

    Lucee::ConstFieldPtr<double> fldPtr = fld.createConstPtr();
    int idx[NDIM];
    double localInt = 0.0;
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      fld.setPtr(fldPtr, idx);
      grid.setIndex(idx);
      double vol = grid.getVolume();
      localInt += vol*getIntegrand(fldPtr);
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
  IntegrateField<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }

  template <unsigned NDIM>
  double
  IntegrateField<NDIM>::getIntegrand(const Lucee::ConstFieldPtr<double>& inp)
  {
    if (!hasFunction)
      return inp[0];

// process field through Lua function
    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;
// push function object on stack
    lua_rawgeti(*L, LUA_REGISTRYINDEX, fnRef);
// push variables on stack
    for (unsigned i=0; i<inp.getNumComponents(); ++i)
      lua_pushnumber(*L, inp[i]);
// call function
    if (lua_pcall(*L, inp.getNumComponents(), 1, 0) != 0)
    {
      std::string err(lua_tostring(*L, -1));
      lua_pop(*L, 1);
      Lucee::Except lce("IntegrateField::getIntegrand: ");
      lce << "Problem evaluating function supplied as 'integrand' ";
      lce << std::endl << "[" << err << "]";
      throw lce;
    }
// fetch results
    if (!lua_isnumber(*L, -1))
        throw Lucee::Except("IntegrateField::getIntegrand: value not a number");
    double res = lua_tonumber(*L, -1);
    lua_pop(*L, 1);

    return res;
  }

// instantiations
  template class Lucee::IntegrateField<1>;
  template class Lucee::IntegrateField<2>;
  template class Lucee::IntegrateField<3>;
}

