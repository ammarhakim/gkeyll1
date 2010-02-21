/**
 * @file	LcFieldPtr.cpp
 *
 * @brief	Pointer to values stored in a field.
 *
 * @version	$Id: LcFieldPtr.cpp 222 2009-11-17 04:46:22Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// lucee includes
#include <LcFieldPtr.h>

namespace Lucee
{
  template <typename T>
  FieldPtr<T>::FieldPtr(unsigned nc, T *data)
    : numComponents(nc), data(data) 
  {
  }

// instantiations
  template class Lucee::FieldPtr<int>;
  template class Lucee::FieldPtr<float>;
  template class Lucee::FieldPtr<double>;
}
