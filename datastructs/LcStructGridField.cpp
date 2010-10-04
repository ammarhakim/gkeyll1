/**
 * @file	LcStructGridField.cpp
 *
 * @brief	StructGridFields are fields that live on structured grids.
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
#include <LcFieldFactory.h>
#include <LcStructGridField.h>

namespace Lucee
{
// names used in registration system
  template <> const char *StructGridField<1, double>::id = "Field1D";
  template <> const char *StructGridField<2, double>::id = "Field2D";
  template <> const char *StructGridField<3, double>::id = "Field3D";
  template <> const char *StructGridField<4, double>::id = "Field4D";
  template <> const char *StructGridField<5, double>::id = "Field5D";
  template <> const char *StructGridField<6, double>::id = "Field6D";
  template <> const char *StructGridField<7, double>::id = "Field7D";

  template <unsigned NDIM, typename T>
  StructGridField<NDIM, T>::StructGridField()
    : Lucee::Field<NDIM, T>(), grid(0)
  {
  }

  template <unsigned NDIM, typename T>
  StructGridField<NDIM, T>::StructGridField(const StructGridField<NDIM, T>& fld)
    : Lucee::Field<NDIM, T>(fld), grid(fld.grid)
  {
  }

  template <unsigned NDIM, typename T>
  StructGridField<NDIM, T>&
  StructGridField<NDIM, T>::operator=(const StructGridField<NDIM, T>& fld)
  {
    if (&fld == this)
      return *this;
// call base class assignment operator
    Field<NDIM, T>::operator=(fld);
    grid = fld.grid;

    return *this;
  }

  template <unsigned NDIM, typename T>
  void
  StructGridField<NDIM, T>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::FieldFactory<NDIM, T> ff;
    ff.readInput(tbl);
// re-rest ourself from factory produced grid
    Field<NDIM, T>::operator=(*ff.create());
    grid = ff.getGridPtr(); // set out grid pointer
  }

  template <unsigned NDIM, typename T>
  Lucee::IoNodeType
  StructGridField<NDIM, T>::writeToFile(Lucee::IoBase& io, Lucee::IoNodeType& node, const std::string& nm)
  {
// first write the field data to file
    Lucee::IoNodeType dn 
      = Lucee::Field<NDIM, T>::writeToFile(io, node, "StructGridField");
// annotate with viz-schema marks
    io.writeStrAttribute(dn, "vsType", "variable");
    io.writeStrAttribute(dn, "vsMesh", "StructGrid");
    io.writeStrAttribute(dn, "vsCentering", "zonal");
// now write out grid
    grid->writeToFile(io, node, "StructGrid");

    return node;
  }

  template <unsigned NDIM, typename T>
  void
  StructGridField<NDIM, T>::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
    lfm.appendFunc("set", luaSet);
  }

  template <unsigned NDIM, typename T>
  int
  StructGridField<NDIM, T>::luaSet(lua_State *L)
  {

    return 0;
  }

// instantiations
  template class StructGridField<1, int>;
  template class StructGridField<2, int>;
  template class StructGridField<3, int>;
  template class StructGridField<4, int>;
  template class StructGridField<5, int>;
  template class StructGridField<6, int>;
  template class StructGridField<7, int>;

  template class StructGridField<1, float>;
  template class StructGridField<2, float>;
  template class StructGridField<3, float>;
  template class StructGridField<4, float>;
  template class StructGridField<5, float>;
  template class StructGridField<6, float>;
  template class StructGridField<7, float>;

  template class StructGridField<1, double>;
  template class StructGridField<2, double>;
  template class StructGridField<3, double>;
  template class StructGridField<4, double>;
  template class StructGridField<5, double>;
  template class StructGridField<6, double>;
  template class StructGridField<7, double>;
}
