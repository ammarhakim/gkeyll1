/**
 * @file	LcDynVector.cpp
 *
 * @brief	Hold a set of global time-dependent values.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDynVector.h>

namespace Lucee
{
  template <> const char *DynVector<double>::id = "DynVector";

  template <typename T>
  DynVector<T>::DynVector()
    : numComponents(1) 
  {
  }

  template <typename T>
  DynVector<T>::DynVector(unsigned nc)
    : numComponents(nc) 
  {
  }

  template <typename T>
  void
  DynVector<T>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    DataStructIfc::readInput(tbl);

    numComponents = 1;
// read in number of components
    if (tbl.hasNumber("numComponents"))
      numComponents = (unsigned) tbl.getNumber("numComponents");
  }

  template <typename T>
  unsigned
  DynVector<T>::getNumComponents() const
  {
    return numComponents;
  }

  template <typename T>
  unsigned
  DynVector<T>::getSize() const
  {
    return data.size();
  }

  template <typename T>
  void
  DynVector<T>::clear() const
  {
    data.clear();
    timeMesh.clear();
  }

  template <typename T>
  void
  DynVector<T>::appendData(double t, const std::vector<T>& dat)
  {
    if (dat.size() != numComponents) {
      Lucee::Except lce("Lucee::DynVector::appendData: Data size should be ");
      lce << numComponents << ". Was " << dat.size() << " instead";
      throw lce;
    }
    lastData = dat;
    data.push_back(dat);
    timeMesh.push_back(t);
  }

  template <typename T>
  std::vector<T>
  DynVector<T>::getLastInsertedData() const
  {
    return lastData;
  }

  template <typename T>
  double
  DynVector<T>::getLastInsertedTime() const
  {
    unsigned n = timeMesh.size();
    if (n>0)
      return timeMesh[n-1];
    throw Lucee::Except("Lucee::DynVector::getLastInsertedTime: Dynamic vector is empty.");
  }

  template <typename T>
  void
  DynVector<T>::removeLastInsertedData()
  {
    data.pop_back();
    timeMesh.pop_back();
  }

  template <typename T>
  TxIoNodeType
  DynVector<T>::writeTimeMesh(const TxIoBase& io,
    const TxIoNodeType& node) const 
  {
// space for passing to HDF I/O wrappers
    std::vector<size_t> dataSetSize;
    std::vector<size_t> dataSetBeg;
    std::vector<size_t> dataSetLen;

// create arrays to pass to HDF5
    dataSetSize.push_back(getSize());
    dataSetBeg.push_back(0);
    dataSetLen.push_back(getSize());

    dataSetSize.push_back(1);
    dataSetBeg.push_back(0);
    dataSetLen.push_back(1);

// write data
    TxIoNodeType dn
      = io.writeDataSet(node, "timeMesh",
        dataSetSize, dataSetBeg, dataSetLen, &timeMesh[0]);

// following attributes are picked up by the Viz tool
    io.writeAttribute(dn, "vsType", "mesh");
    io.writeAttribute(dn, "vsKind", "structured");

    return dn;
  }

  template <typename T>
  TxIoNodeType
  DynVector<T>::writeData(const TxIoBase& io,
    const TxIoNodeType& node) const 
  {
// space for passing to HDF I/O wrappers
    std::vector<size_t> dataSetSize;
    std::vector<size_t> dataSetBeg;
    std::vector<size_t> dataSetLen;

// create arrays to pass to HDF5
    unsigned vol = getNumComponents()*getSize();
    std::vector<T> rdata(vol);
    dataSetSize.push_back(getSize());
    dataSetBeg.push_back(0);
    dataSetLen.push_back(getSize());

    dataSetSize.push_back(getNumComponents());
    dataSetBeg.push_back(0);
    dataSetLen.push_back(getNumComponents());

// copy data into array
    unsigned count = 0;
    typename std::vector<std::vector<T> >::const_iterator itr;
    for (itr = data.begin(); itr != data.end(); ++itr) {
      for (unsigned n=0; n<getNumComponents(); ++n) {
        rdata[count++] = itr->at(n);
      }
    }
// write it
    TxIoNodeType dn = io.writeDataSet(node, "data",
      dataSetSize, dataSetBeg, dataSetLen, &rdata[0]);

// following attributes are picked up by the Viz tool
    io.writeAttribute(dn, "vsType", "variable");
    io.writeAttribute(dn, "vsMesh", "timeMesh");

    return dn;
  }

  template <typename T>
  TxIoNodeType
  DynVector<T>::writeToFile(TxIoBase& io, TxIoNodeType& node,
    const std::string& name)
  {
// skip writing if nothing to write
    if (data.size() == 0)
      return node;

// create a new group to write data
    TxIoNodeType grp = io.createGroup(node, name);

// now write time-mesh and data to this group
    writeTimeMesh(io, grp);
    writeData(io, grp);

// now clear all data before next write
    clear();

    return node;
  }

  template <typename T>
  void
  DynVector<T>::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
    DataStructIfc::appendLuaCallableMethods(lfm);
    lfm.appendFunc("lastInsertedData", luaLastInsertedData);
    lfm.appendFunc("lastInsertedTime", luaLastInsertedTime);
  }

  template <typename T>
  int
  DynVector<T>::luaLastInsertedData(lua_State *L)
  {
    DynVector<T> *dynv
      = Lucee::PointerHolder<DynVector<T> >::getObjAsDerived(L);

    std::vector<T> data = dynv->getLastInsertedData();
    for (unsigned i=0; i<data.size(); ++i)
      lua_pushnumber(L, data[i]);
    return dynv->getNumComponents();
  }

  template <typename T>
  int
  DynVector<T>::luaLastInsertedTime(lua_State *L)
  {
    DynVector<T> *dynv
      = Lucee::PointerHolder<DynVector<T> >::getObjAsDerived(L);
    lua_pushnumber(L, dynv->getLastInsertedTime());
    return 1;
  }

// instantiations
  template class DynVector<float>;
  template class DynVector<double>;
}
