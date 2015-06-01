/**
 * @file	LcUnstructGridField.h
 *
 * @brief	UnstructGridFields are fields that live on structured grids.
 */

#ifndef LC_UNSTRUCT_GRID_FIELD_H
#define LC_UNSTRUCT_GRID_FIELD_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcQuadratureRule.h>
#include <LcStructuredGridBase.h>

// txbase includes
#include <TxCommBase.h>

namespace Lucee
{
/**
 * A structured field is a field that lives on a structured grid.
 */
  template <unsigned NDIM, typename T>
  class UnstructGridField : public Lucee::Field<NDIM, T>
  {
    public:
/** Class id: this is used by the registration system */
      static const char *id;

/**
 * Create an empty field. This should not be use directly: it is
 * provided for use in creation of fields from Lua scripts.
 */      
      UnstructGridField();

/**
 * Create a new field indexing given region. This constructor creates
 * an empty set of ghost indices.
 *
 * @param grid Global grid on which field lives.
 * @param nc Number of components at each index location.
 * @param lg Ghost indexes along lower index range in each dimension.
 * @param ug Ghost indexes along upper index range in each dimension.
 */      
      UnstructGridField(Lucee::StructuredGridBase<NDIM>* grid, unsigned nc,
        int lg[NDIM], int ug[NDIM]);

/**
 * Create a field from supplied one (shallow copy).
 *
 * @param fld UnstructGridField to copy from.
 */
      UnstructGridField(const UnstructGridField<NDIM, T>& fld);

/**
 * Assign a field from supplied one (shallow copy).
 *
 * @param fld UnstructGridField to assign from.
 * @return Reference to this field.
 */
      UnstructGridField<NDIM, T>& operator=(const UnstructGridField<NDIM, T>& fld);

/**
 * Set all values of field to supplied one.
 *
 * @param val Value to set.
 * @return Reference to this field.
 */
      UnstructGridField<NDIM, T>& operator=(const T& val);

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
      void divergence(Lucee::UnstructGridField<NDIM, T>& div) const;

/**
 * Write dataStruct to given node in HDF5 file.
 *
 * @param io I/O object for I/O.
 * @param node Node to write to.
 * @param nm Name of the array as it should appear in output.
 * @return node to which data was written.
 */
      virtual TxIoNodeType writeToFile(TxIoBase& io, TxIoNodeType& node,
         const std::string& nm);

/**
 * Read dataStruct from given node in HDF5 file.
 *
 * @param io I/O object for I/O.
 * @param node Node to read data from.
 * @param nm Name of the data-struct as it appears in input
 * @return node from which data was read.
 */
      virtual TxIoNodeType readFromFile(TxIoBase& io, TxIoNodeType& node,
        const std::string& nm);

/**
 * Write dataStruct to specified text file. The data is written as
 * plain text.
 *
 * @param txtFl Text file handle for output.
 */
      virtual void writeToTxtFile(std::ofstream& txtFl);

/**
 * Synchronise data in ghost cells, copying skin data into neighbor
 * ghost cells.
 */
      virtual void sync();

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
 * Lua callable method for initializing ghost region of field from Lua
 * supplied function.
 *
 * @param L Lua state to use.
 * @return number of output parameters.
 */
      static int luaSetGhost(lua_State *L);

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

/**
 * Set values in the ghost region of the field from Lua function. The
 * function itself is specified using a reference to the Lua function
 * on stack.
 *
 * @param L Lua state to use.
 * @param ref Reference to Lua function.
 * @param dir Direction to apply function.
 * @param side Side to apply function (0 - lower, 1 - upper).
 */
      void setGhostFromLuaFunction(lua_State *L, int ref, unsigned dir, unsigned side);

    private:
/** Pointer to grid */
      Lucee::StructuredGridBase<NDIM> *grid;
/** Flag to indicate location of data in field */
      int dataLoc;
/** Flag to indicate if we are in recieve mode */
      mutable bool isReceiving;
/** Map of a rank -> pending recieves */
      mutable std::map<int, TxMsgStatus> msgStatus;
/** Lower ghost cells to write */
      int lowerWriteGhost[NDIM];
/** Upper ghost cells to write */
      int upperWriteGhost[NDIM];

/**
 * Start a non-blocking operation to receive data from skin region.
 *
 * @param lowerExt Length of extension along lower end in each direction.
 * @param upperExt Length of extension along upper end in each direction.
 */
      void startRecv(unsigned rank, int lowerGhost[NDIM], int upperGhost[NDIM]);

/**
 * Complete the non-blocking operation to receive data from skin region.
 */
    void finishRecv();

/**
 * Do a blocking operation to send data to our neighbors.
 *
 * @param lowerExt Length of extension along lower end in each direction.
 * @param upperExt Length of extension along upper end in each direction.
 */
    void send(unsigned rank, int lowerGhost[NDIM], int upperGhost[NDIM]);

/**
 * Create a field from supplied one (shallow copy).
 *
 * @param fld Field to copy from.
 * @param grd Grid to use
 */
      UnstructGridField(const Field<NDIM, T>& fld, Lucee::StructuredGridBase<NDIM>& grd);
  };
}

#endif // LC_STRUCT_GRID_FIELD_H
