/**
 * @file	LcLuaObjTypeId.cpp
 *
 * @brief       Class that identifies an object derived from BasicObj.
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
  LuaObjTypeId::checkDerivedTypeId(lua_State *L, const std::string& dtype, void **obj, int loc)
  {
    Lucee::PointerHolder<Lucee::BasicObj> *ph =
      (Lucee::PointerHolder<Lucee::BasicObj>*) lua_touserdata(L, loc);

    if (ph->pointer->getDerivedType() == dtype)
    {
      *obj = (void*) ph;
      return true;
    }
// did not match
    *obj = 0;
    return false;
  }

  bool
  LuaObjTypeId::checkBaseTypeId(lua_State *L, const std::string& btype, void **obj, int loc)
  {
    Lucee::PointerHolder<Lucee::BasicObj> *ph =
      (Lucee::PointerHolder<Lucee::BasicObj>*) lua_touserdata(L, loc);

    if (ph->pointer->getBaseType() == btype)
    {
      *obj = (void*) ph;
      return true;
    }
// did not match
    *obj = 0;
    return false;
  }
}
