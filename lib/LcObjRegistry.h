/**
 * @file	LcObjRegistry.h
 *
 * @brief	Class for handling object registration
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_OBJ_REGISTRY_H
#define LC_OBJ_REGISTRY_H

// lucee includes
#include <LcExcept.h>
#include <LcLuaModule.h>
#include <LcLuaTable.h>
#include <LcPointerHolder.h>

// loki includes
#include <loki/Factory.h>
#include <loki/Singleton.h>

// std includes
#include <typeinfo>

namespace Lucee
{
/**
 * Class to perform object registration.
 */
  template<class B, class D>
  class ObjRegistry
  {
    public:
/**
 * Register a new object to be created by its id.
 */
      ObjRegistry()
      {
        Lucee::LuaModule<B>& lm = Loki::SingletonHolder<Lucee::LuaModule<B> >
          ::Instance();
// add a function to make Lua object from table constructor
        luaL_Reg reg = {D::id, makeLuaObj};
        lm.regCreateFuncs.push_back(reg);

// add an empty map to using derived class ID as key
        lm.funcMaps.insert(std::pair<std::string, Lucee::LuaFuncMap>(
            typeid(D).name(), Lucee::LuaFuncMap()));
        std::map<std::string, Lucee::LuaFuncMap>::iterator itr
          = lm.funcMaps.find(typeid(D).name()); // fetch just-added map entry
// set its base class name
        itr->second.setBaseName(typeid(B).name());
// add Lua callable functions for base class
        B::appendLuaCallableMethods(itr->second);
// add Lua callable functions for derived class
        D::appendLuaCallableMethods(itr->second);
// set deletion function for derived class
        itr->second.setDelFunc(Lucee::PointerHolder<D>::deleteObject);
      }

/**
 * Function called by Lua for constructing object and returning it as
 * a Lua variable.
 *
 * @return Newly allocated object.
 */
      static int makeLuaObj(lua_State *L)
      {
        size_t nbytes = sizeof(Lucee::PointerHolder<D>);
        Lucee::PointerHolder<D> *ph =
          (Lucee::PointerHolder<D>*) lua_newuserdata(L, nbytes);
// make the table top of stack
        lua_pushvalue(L, 1);
        Lucee::LuaState myL(L);
        Lucee::LuaTable tbl(myL, "Module");
        lua_pop(L, 1); // pop off table from top

// create object
        ph->pointer = new D;
// initialize object using Lua table
        ph->pointer->readInput(tbl);
// set base and derived types
        ph->pointer->template setBaseType<B>();
        ph->pointer->template setDerivedType<D>();

// get a meta-table for this object and set it
        luaL_getmetatable(L, typeid(D).name());
        lua_setmetatable(L, -2);

        return 1;
      }
  };
}

#endif // LC_OBJ_REGISTRY_H
