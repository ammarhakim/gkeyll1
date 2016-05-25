/**
 * @file	LcFieldFunctionUpdater.cpp
 *
 * @brief	Field function updater / 'calculator'
 */

// lucee includes
#include <LcFieldFunctionUpdater.h>
#include <LcStructGridField.h>
#include <LcGlobals.h>
#include <LcLuaState.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{

// set ids for module system
  template <> const char *FieldFunctionUpdater<1>::id = "FieldFunction1D";
  template <> const char *FieldFunctionUpdater<2>::id = "FieldFunction2D";
  template <> const char *FieldFunctionUpdater<3>::id = "FieldFunction3D";

  template <unsigned NDIM>
  void
  FieldFunctionUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

// get list of components to apply this boundary condition
    std::vector<double> cd = tbl.getNumVec("inpComponents");
    for (unsigned i=0; i<cd.size(); ++i)
      inpComponents.push_back( (int) cd[i] );

// get list of components to apply this boundary condition
    cd = tbl.getNumVec("outComponents");
    for (unsigned i=0; i<cd.size(); ++i)
      outComponents.push_back( (int) cd[i] );

// get reference to function
    fnRef = tbl.getFunctionRef("func");

  }

  template <unsigned NDIM>
  void
  FieldFunctionUpdater<NDIM>::initialize()
  {
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  FieldFunctionUpdater<NDIM>::update(double t)
  {
    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;

    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    const Lucee::Field<NDIM, double>& inpField = this->getInp<Lucee::Field<NDIM, double> >(0);
    Lucee::Field<NDIM, double>& outField = this->getOut<Lucee::Field<NDIM, double> >(0);

    Lucee::ConstFieldPtr<double> inpPtr = inpField.createConstPtr();
    Lucee::FieldPtr<double> outPtr = outField.createPtr();

    int idx[NDIM];
    double xc[3];
    Lucee::Region<NDIM, int> localRgn = inpField.getRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      
      grid.setIndex(idx);
      grid.getCentroid(xc);

      inpField.setPtr(inpPtr, idx);
      outField.setPtr(outPtr, idx);

// push function object on stack
      lua_rawgeti(*L, LUA_REGISTRYINDEX, fnRef);
// push variables on stack
      for (unsigned i=0; i<3; ++i)
        lua_pushnumber(*L, xc[i]);
      lua_pushnumber(*L, t);
// push
      for (unsigned i=0; i<inpComponents.size(); ++i)
        lua_pushnumber(*L, inpPtr[inpComponents[i]]);

      unsigned numInp = 4+inpComponents.size();
      unsigned numOut = outComponents.size();
// call function
      if (lua_pcall(*L, numInp, numOut, 0) != 0)
      {
        std::string err(lua_tostring(*L, -1));
        lua_pop(*L, 1);
        Lucee::Except lce("FieldFunctionUpdater::update: ");
        lce << "Problem evaluating function supplied as 'func' ";
        lce << std::endl << "[" << err << "]";
        throw lce;
      }
// fetch results
      for (int i=-outComponents.size(); i<0; ++i)
      {
        if (!lua_isnumber(*L, i))
          throw Lucee::Except("FieldFunctionUpdater::update: Return value not a number");
        outPtr[outComponents[outComponents.size()+i]] = lua_tonumber(*L, i);
      }
      lua_pop(*L, 1);
    }
    
    return Lucee::UpdaterStatus();
  }
  

  template <unsigned NDIM>
  void
  FieldFunctionUpdater<NDIM>::declareTypes()
  {
    this->setLastOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class FieldFunctionUpdater<1>;
  template class FieldFunctionUpdater<2>;
  template class FieldFunctionUpdater<3>;
}

