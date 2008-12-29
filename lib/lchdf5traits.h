/**
 * @file	lchdf5traits.h
 *
 * @brief	Class to define traits for types supported by HDF5.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_HDF5_TRAITS_H
#define LC_HDF5_TRAITS_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

/** Hdf5 */
#define H5_USE_16_API
#include <hdf5.h>

namespace Lucee
{
/**
 * Traits class for HDF5
 */
  template <typename T> struct Hdf5Traits;

  template <> struct Hdf5Traits<char> 
  {
/**
 * Return HDF5 type for this C++ type.
 * 
 * @return HDF5 type.
 */
      static hid_t hdf5Type() 
      {
        return H5T_NATIVE_CHAR;
      }
  };

  template <> struct Hdf5Traits<unsigned char> 
  {
/**
 * Return HDF5 type for this C++ type.
 * 
 * @return HDF5 type.
 */
      static hid_t hdf5Type() 
      {
        return H5T_NATIVE_UCHAR;
      }
  };
  
  template <> struct Hdf5Traits<short> 
  {
/**
 * Return HDF5 type for this C++ type.
 * 
 * @return HDF5 type.
 */
      static hid_t hdf5Type() 
      {
        return H5T_NATIVE_SHORT;
      }
  };

  template <> struct Hdf5Traits<unsigned short> 
  {
/**
 * Return HDF5 type for this C++ type.
 * 
 * @return HDF5 type.
 */
      static hid_t hdf5Type() 
      {
        return H5T_NATIVE_USHORT;
      }
  };

  template <> struct Hdf5Traits<int> 
  {
/**
 * Return HDF5 type for this C++ type.
 * 
 * @return HDF5 type.
 */
      static hid_t hdf5Type() 
      {
        return H5T_NATIVE_INT;
      }
  };

  template <> struct Hdf5Traits<unsigned int> 
  {
/**
 * Return HDF5 type for this C++ type.
 * 
 * @return HDF5 type.
 */
      static hid_t hdf5Type() 
      {
        return H5T_NATIVE_UINT;
      }
  };

  template <> struct Hdf5Traits<long> 
  {
/**
 * Return HDF5 type for this C++ type.
 * 
 * @return HDF5 type.
 */
      static hid_t hdf5Type() 
      {
        return H5T_NATIVE_LONG;
      }
  };

  template <> struct Hdf5Traits<unsigned long> 
  {
/**
 * Return HDF5 type for this C++ type.
 * 
 * @return HDF5 type.
 */
      static hid_t hdf5Type() 
      {
        return H5T_NATIVE_ULONG;
      }
  };

  template <> struct Hdf5Traits<float> 
  {
/**
 * Return HDF5 type for this C++ type.
 * 
 * @return HDF5 type.
 */
      static hid_t hdf5Type() 
      {
        return H5T_NATIVE_FLOAT;
      }
  };

  template <> struct Hdf5Traits<double> 
  {
/**
 * Return HDF5 type for this C++ type.
 * 
 * @return HDF5 type.
 */
      static hid_t hdf5Type() 
      {
        return H5T_NATIVE_DOUBLE;
      }
  };

  template <> struct Hdf5Traits<long double>  
  {
/**
 * Return HDF5 type for this C++ type.
 * 
 * @return HDF5 type.
 */
      static hid_t hdf5Type() 
      {
        return H5T_NATIVE_LDOUBLE;
      }
  };

  template <> struct Hdf5Traits<long long> 
  {
/**
 * Return HDF5 type for this C++ type.
 * 
 * @return HDF5 type.
 */
      static hid_t hdf5Type() 
      {
        return H5T_NATIVE_LLONG;
      }
  };
}

#endif // LC_HDF5_TRAITS_H
