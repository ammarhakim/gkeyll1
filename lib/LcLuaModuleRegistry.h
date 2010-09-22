/**
 * @file	LcLuaModuleRegistry.h
 *
 * @brief	Class for handling object creation.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_LUA_MODULE_REGISTRY_H
#define LC_LUA_MODULE_REGISTRY_H

// lucee includes
#include <LcExcept.h>
#include <LcLuaModule.h>

// loki includes
#include <loki/Factory.h>
#include <loki/Singleton.h>

namespace Lucee
{
/**
 * Class to create objects. This class is used in conjunction with the
 * Lucee::ObjRegistry class. Once an object has been registered this
 * class can be used to create derived classes using their names.
 */
  template <class B>
  class LuaModuleRegistry
  {
    public:
/**
 * Register modules. This function will register Lua functions to
 * construct Lucee objects. The functions will become avaiable for use
 * from Lua scripts.
 *
 * @param L Lua state in which modules should be registered.
 */
      static void registerModule(lua_State *L)
      {
// create name for a metatable for this object type
        std::string mtblNm("Lucee.");
        mtblNm.append(std::string(B::id));
        luaL_newmetatable(L, mtblNm.c_str()); // create a new metatable
// register the functions to create children objects
        luaL_Reg reg = {NULL, NULL};
        Loki::SingletonHolder<Lucee::LuaModule<B> >
          ::Instance().regObjFuncs.push_back(reg);
        luaL_register(L, B::id, &Loki::SingletonHolder<Lucee::LuaModule<B> >
          ::Instance().regObjFuncs[0]);
      }

/**
 * Get a list of registered names.
 *
 * @return List of registered names.
 */
      static std::vector<std::string> registeredNames() 
      {
        return Loki::SingletonHolder<Loki::Factory<B, std::string> >
          ::Instance().RegisteredIds();
      }

/**
 * Unregister creator with given name.
 *
 * @param nm Name of creator to unregister.
 * @return bool true if unregister worked, false otherwise.
 */
      static bool unregister(const std::string& nm)
      {
        return Loki::SingletonHolder<Loki::Factory<B, std::string> >
          ::Instance().Unregister(nm);
      }
  };
}

#endif // LC_LUA_MODULE_REGISTRY_H
