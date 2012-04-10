/**
 * @file	LcDynVector.cpp
 *
 * @brief	Hold a set of global time-depent values.
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
  void
  DynVector<T>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    DataStructIfc::readInput(tbl);
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
  DynVector<T>::writeToFile(TxIoBase& io, TxIoNodeType& node,
    const std::string& nm)
  {
  }

// instantiations
  template class DynVector<double>;
}
