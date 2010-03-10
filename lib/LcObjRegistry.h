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
// add a function
        luaL_Reg reg = {D::id, getModuleTable};
        Loki::SingletonHolder<Lucee::LuaModule<B> >
          ::Instance().regFuncs.push_back(reg);
      }

/**
 * Function called by Lua for constructing object table.
 *
 * @return Newly allocated object.
 */
      static int getModuleTable(lua_State *L)
      {
        if (! lua_istable(L, 1) )
          throw Lucee::Except("ObjRegistry::getModuleTable: expects table as single parameter");
// add meta-data keys into table
        lua_pushstring(L, "__type");
        lua_pushstring(L, B::id);
        lua_settable(L, -3);

        lua_pushstring(L, "__kind");
        lua_pushstring(L, D::id);
        lua_settable(L, -3);

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
