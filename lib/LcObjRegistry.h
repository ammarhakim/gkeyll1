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
// register creator function
        Loki::SingletonHolder<Loki::Factory<B, std::string> >
          ::Instance().Register(D::id, getNew);
// add a function to make Lua object
        luaL_Reg reg = {D::id, makeLuaObj};
        Loki::SingletonHolder<Lucee::LuaModule<B> >
          ::Instance().regObjFuncs.push_back(reg);
      }

/**
 * Function called by Lua for constructing object table.
 *
 * @return Newly allocated object.
 */
      static int makeLuaObj(lua_State *L)
      {
        size_t nbytes = sizeof(Lucee::PointerHolder<B>);
        Lucee::PointerHolder<B> *ph =
          (Lucee::PointerHolder<B>*) lua_newuserdata(L, nbytes);
// make the table top of stack
        lua_pushvalue(L, 1);
        Lucee::LuaState myL(L);
        Lucee::LuaTable tbl(myL, "Module");
        lua_pop(L, 1); // pop off table from top

// create object
        ph->pointer = getNew();
// setup object
        ph->pointer->readInput(tbl);

        return 1;
      }

    private:
/**
 * Return a newly allocated object. This is used as the creation
 * mechanism for the object.
 *
 * @return Newly allocated object.
 */
      static B* getNew() 
      {
        return new D;
      }
  };
}

#endif // LC_OBJ_REGISTRY_H
