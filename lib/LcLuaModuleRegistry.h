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
#include <LcPointerHolder.h>

// loki includes
#include <loki/Factory.h>
#include <loki/Singleton.h>

namespace Lucee
{
/**
 * Class to register Lua callable functions. The registerModule() method
 * must be called for each base class type to register the table
 * constructors for all its derived classes. In addition, this method
 * also registers the Lua-callable functions
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
// register the functions to create children objects
        luaL_Reg reg = {NULL, NULL};
        Loki::SingletonHolder<Lucee::LuaModule<B> >
          ::Instance().regCreateFuncs.push_back(reg);

        Lucee::LuaModule<B>& lm = Loki::SingletonHolder<Lucee::LuaModule<B> >
          ::Instance();
        luaL_register(L, B::id, &lm.regCreateFuncs[0]);

// now create a meta-table for each derived class methods and register
// them so that they become available using Lua OO notation.
        std::map<std::string, Lucee::LuaFuncMap>::iterator itr
          = lm.funcMaps.begin();
        for ( ; itr != lm.funcMaps.end(); ++itr)
        {
          luaL_newmetatable(L, itr->first.c_str());
          lua_pushvalue(L, -1); // copy metatable
          lua_setfield(L, -2, "__index");
// get list of functions
          std::vector<luaL_Reg> funcLst;
          itr->second.fillWithFuncList(funcLst);
// append delete so object is cleaned-up when garbage collector runs
          luaL_Reg gc = {"__gc", itr->second.getDelFunc()};
          funcLst.push_back(gc);
          luaL_Reg fin = {NULL, NULL};
          funcLst.push_back(fin);
// now register it
          luaL_register(L, NULL, &funcLst[0]);
        }
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
