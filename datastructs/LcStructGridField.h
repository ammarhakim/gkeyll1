/**
 * @file	LcStructGridField.h
 *
 * @brief	StructGridFields are fields that live on structured grids.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_STRUCT_GRID_FIELD_H
#define LC_STRUCT_GRID_FIELD_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcStructuredGridBase.h>

namespace Lucee
{
/**
 * A structured field is a field that lives on a structured grid.
 */
  template <unsigned NDIM, typename T>
  class StructGridField : public Lucee::Field<NDIM, T>
  {
    public:
/** Class id: this is used by the registration system */
      static const char *id;

/**
 * Create an empty field. This should not be use directly: it is
 * provided for use in creation of fields from Lua scripts.
 */      
      StructGridField();

/**
 * Create a field from supplied one (shallow copy).
 *
 * @param fld StructGridField to copy from.
 */
      StructGridField(const StructGridField<NDIM, T>& fld);

/**
 * Assign a field from supplied one (shallow copy).
 *
 * @param fld StructGridField to assign from.
 * @return Reference to this field.
 */
      StructGridField<NDIM, T>& operator=(const StructGridField<NDIM, T>& fld);

/**
 * Set all values of field to supplied one.
 *
 * @param val Value to set.
 * @return Reference to this field.
 */
      StructGridField<NDIM, T>& operator=(const T& val);

/**
 * Create from Lua table data.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Compute divergence of this field and store in supplied field.
 *
 * @param div Divergence is stored in this field.
 */
      void divergence(Lucee::StructGridField<NDIM, T>& div) const;

/**
 * Write dataStruct to given node in HDF5 file.
 *
 * @param io I/O object for I/O.
 * @param node Node to write to.
 * @param nm Name of the array as it should appear in output.
 * @return node to which data was written.
 */
      virtual Lucee::IoNodeType writeToFile(Lucee::IoBase& io, Lucee::IoNodeType& node,
         const std::string& nm);

/**
 * Write dataStruct to specified text file. The data is written as
 * plain text.
 *
 * @param txtFl Text file handle for output.
 */
      virtual void writeToTxtFile(std::ofstream& txtFl);

/**
 * Method that performs registration of Lua functions.
 *
 * @param lfm Lua function map object.
 */
      static void appendLuaCallableMethods(Lucee::LuaFuncMap& lfm);

/**
 * Lua callable method for initializing field from Lua supplied
 * function.
 *
 * @param L Lua state to use.
 * @return number of output parameters.
 */
      static int luaSet(lua_State *L);

/**
 * Lua callable method to return an alias to this field.
 *
 * @param L Lua state to use.
 * @return number of output parameters.
 */
      static int luaAlias(lua_State *L);

/**
 * Lua callable method to return duplicate of this field.
 *
 * @param L Lua state to use.
 * @return number of output parameters.
 */
      static int luaDuplicate(lua_State *L);

/**
 * Lua callable method to compute divergence of this field.
 *
 * @param L Lua state to use.
 * @return number of output parameters.
 */
      static int luaDivergence(lua_State *L);

/**
 * Set field from Lua function. The function itself is specified using
 * a reference to the Lua function on stack.
 *
 * @param L Lua state to use.
 * @param ref Reference to Lua function.
 */
      void setFromLuaFunction(lua_State *L, int ref);

    private:
/** Pointer to grid */
      Lucee::StructuredGridBase<NDIM> *grid;

/**
 * Create a field from supplied one (shallow copy).
 *
 * @param fld Field to copy from.
 * @param grd Grid to use
 */
      StructGridField(const Field<NDIM, T>& fld, Lucee::StructuredGridBase<NDIM>& grd);
  };
}

#endif // LC_STRUCT_GRID_FIELD_H
