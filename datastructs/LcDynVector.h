/**
 * @file	LcDynVector.h
 *
 * @brief	Hold a set of global time-depent values.
 */

#ifndef LC_TIME_VECTOR_H
#define LC_TIME_VECTOR_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDataStructIfc.h>

namespace Lucee
{
/**
 * Datastructure to hold a set of global values that depend on time.
 */
  template <typename T>
  class DynVector : public Lucee::DataStructIfc
  {
    public:
/** Class id: this is used by the registration system */
      static const char *id;

/**
 * Create a new vector having supplied number of components
 *
 * @param numComponents number of components in array
 */
      DynVector();

/**
 * Create a new vector having supplied number of components
 *
 * @param numComponents number of components in array
 */
      DynVector(unsigned numComponents);

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Return number of components per entry in vector.
 *
 * @return components in vector
 */
    unsigned getNumComponents() const;

/**
 * Return number of elements in vector.
 *
 * @return components in vector
 */
    unsigned getSize() const;

/**
 * Clear the data held by this vector.
 */
    void clear() const;

/**
 * Append data to end of this vector.
 *
 * @param t Time at which data was taken.
 * @param data Data vector.
 */
    void appendData(double t, const std::vector<T>& data);

/**
 * Get last data stored in this vector.
 *
 * @return Data vector.
 */
    std::vector<T> getLastInsertedData() const;

/**
 * Get time at which last data value was stored.
 *
 * @return time at which last data value was stored.
 */
    double getLastInsertedTime() const;

/**
 * Remove last data time stored in this vector.
 */
    void removeLastInsertedData();

/**
 * Write dataStruct to given node in HDF5 file.
 *
 * @param io I/O object for I/O.
 * @param node Node to write to.
 * @param nm Name of the data-struct as it should appear in output.
 * @return node to which data was written.
 */
      TxIoNodeType writeToFile(TxIoBase& io, TxIoNodeType& node,
        const std::string& nm);

/**
 * Method that performs registration of Lua functions.
 *
 * @param lfm Lua function map object.
 */
      static void appendLuaCallableMethods(Lucee::LuaFuncMap& lfm);

    private:
/** Number of components to store */
      unsigned numComponents;
/** Data in dyn-vector */
      mutable std::vector<std::vector<T> > data;
/** Time at which data were written */
      mutable std::vector<double> timeMesh;
/** Last stored data */
      std::vector<T> lastData;

/**
 * Write the time-mesh out to h5 file.
 *
 * @param io I/O object for file.
 * @param node Node to write data to.
 * @return Node to which data was written.
 */
    TxIoNodeType writeTimeMesh(const TxIoBase& io, const TxIoNodeType& node) const;

/**
 * Write the data out to h5 file.
 *
 * @param io I/O object for file.
 * @param node Node to write data to.
 * @return Node to which data was written.
 */
    TxIoNodeType writeData(const TxIoBase& io, const TxIoNodeType& node) const;

/**
 * Lua callable method for getting last inserted data.
 *
 * @param L Lua state to use.
 * @return number of output parameters.
 */
      static int luaLastInsertedData(lua_State *L);

/**
 * Lua callable method for getting last inserted time.
 *
 * @param L Lua state to use.
 * @return number of output parameters.
 */
      static int luaLastInsertedTime(lua_State *L);
  };
}

#endif // LC_TIME_VECTOR_H
