/**
 * @file	LcLuaObjTypeId.h
 *
 * @brief       Class that identifies an object derived from BasicObj.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_LUA_OBJ_TYPE_ID_H
#define LC_LUA_OBJ_TYPE_ID_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lua includes
#include <lua.hpp>

// std includes
#include <string>

namespace Lucee
{
/** 
 * Class that allows quering of identifing strings for classes
 * derived from BasicObj. 
 */
  class LuaObjTypeId
  {
    public:
/**
 * Check the type of Lua userdata object against the derived type
 * type-ID 'dtype'. If these match, returns true and sets object
 * pointer as 'obj'. Returns false otherwise and sets obj to 0.
 *
 * @param L Lua state pointer.
 * @param dtype Typeid for derived class.
 * @param obj On output pointer to usedata object.
 */
      static bool checkDerivedTypeId(lua_State *L, const std::string& dtype, void *obj);

/**
 * Check the type of Lua userdata object against the base type type-ID
 * 'dtype'. If these match, returns true and sets object pointer as
 * 'obj'. Returns false otherwise and sets obj to 0.
 *
 * @param L Lua state pointer.
 * @param btype Typeid for base class.
 * @param obj On output pointer to usedata object.
 */
      static bool checkBaseTypeId(lua_State *L, const std::string& btype, void *obj);
  };
}

#endif // LC_LUA_OBJ_TYPE_ID_H
