/**
 * @file	LcFieldPtr.cpp
 *
 * @brief	Pointer to values stored in a field.
 */

// lucee includes
#include <LcFieldPtr.h>

namespace Lucee
{
  template <typename T>
  FieldPtr<T>::FieldPtr(unsigned nc, T *data)
    : numComponents(nc), data(data), isAlloc(false)
  {
  }

  template <typename T>
  FieldPtr<T>::FieldPtr(const FieldPtr<T>& ptr)
    : numComponents(ptr.numComponents), data(ptr.data), isAlloc(ptr.isAlloc)
  {
  }

  template <typename T>
  FieldPtr<T>::FieldPtr(std::vector<T>& vec)
    : numComponents(vec.size()), data(&vec[0]), isAlloc(false)
  {
  }

  template <typename T>
  FieldPtr<T>::FieldPtr(Lucee::Vector<T> vec)
    : numComponents(vec.size()), data(&vec[0]), isAlloc(false)
  {
  }

  template <typename T>
  FieldPtr<T>::FieldPtr(unsigned num)
    : numComponents(num), data(new T[num]), isAlloc(true)
  {
  }

  template <typename T>
  FieldPtr<T>::~FieldPtr()
  {
    if (isAlloc)
      delete [] data;
  }

  template <typename T>
  FieldPtr<T>&
  FieldPtr<T>::operator=(const FieldPtr<T>& ptr)
  {
    if (this == &ptr)
      return *this;

    numComponents = ptr.numComponents;
    data = ptr.data;
    isAlloc = false;
    return *this;
  }

  template <typename T>
  FieldPtr<T>&
  FieldPtr<T>::operator=(const T& val)
  {
    for (unsigned i=0; i<numComponents; ++i)
      data[i] = val;
    return *this;
  }

// instantiations
  template class Lucee::FieldPtr<int>;
  template class Lucee::FieldPtr<float>;
  template class Lucee::FieldPtr<double>;
}
