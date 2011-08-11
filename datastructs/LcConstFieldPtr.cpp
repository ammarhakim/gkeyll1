/**
 * @file	LcConstFieldPtr.cpp
 *
 * @brief	Pointer to values stored in a field.
 */

// lucee includes
#include <LcConstFieldPtr.h>

namespace Lucee
{
  template <typename T>
  ConstFieldPtr<T>::ConstFieldPtr(unsigned nc, const T *data)
    : numComponents(nc), data(data) 
  {
  }

  template <typename T>
  ConstFieldPtr<T>::ConstFieldPtr(const ConstFieldPtr<T>& ptr)
    : numComponents(ptr.numComponents), data(ptr.data)
  {
  }

  template <typename T>
  ConstFieldPtr<T>::ConstFieldPtr(const std::vector<T>& vec)
    : numComponents(vec.size()), data(&vec[0])
  {
  }

  template <typename T>
  ConstFieldPtr<T>::ConstFieldPtr(Lucee::Vector<T> vec)
    : numComponents(vec.size()), data(&vec[0])
  {
  }

  template <typename T>
  ConstFieldPtr<T>::ConstFieldPtr(const FieldPtr<T>& ptr)
    : numComponents(ptr.numComponents), data(ptr.data)
  {
  }

  template <typename T>
  ConstFieldPtr<T>&
  ConstFieldPtr<T>::operator=(const ConstFieldPtr<T>& ptr)
  {
    if (this == &ptr)
      return *this;

    numComponents = ptr.numComponents;
    data = ptr.data;
    return *this;
  }

// instantiations
  template class Lucee::ConstFieldPtr<int>;
  template class Lucee::ConstFieldPtr<float>;
  template class Lucee::ConstFieldPtr<double>;
}
