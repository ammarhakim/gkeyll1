/**
 * @file	LcLuaObjTypeId.cpp
 *
 * @brief       Class that identifies an object derived from BasicObj.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */
// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBasicObj.h>
#include <LcLuaObjTypeId.h>
#include <LcPointerHolder.h>

namespace Lucee
{
  bool
  LuaObjTypeId::checkDerivedTypeId(lua_State *L, const std::string& dtype, void *obj)
  {
    Lucee::PointerHolder<Lucee::BasicObj> *ph =
      (Lucee::PointerHolder<Lucee::BasicObj>*) lua_touserdata(L, 1);

    if (ph->pointer->getDerivedType() == dtype)
    {
      obj = (void*) ph;
      return true;
    }
// did not match
    obj = 0;
    return false;
  }

  bool
  LuaObjTypeId::checkBaseTypeId(lua_State *L, const std::string& btype, void *obj)
  {
    Lucee::PointerHolder<Lucee::BasicObj> *ph =
      (Lucee::PointerHolder<Lucee::BasicObj>*) lua_touserdata(L, 1);

    if (ph->pointer->getBaseType() == btype)
    {
      obj = (void*) ph;
      return true;
    }
// did not match
    obj = 0;
    return false;
  }
}
