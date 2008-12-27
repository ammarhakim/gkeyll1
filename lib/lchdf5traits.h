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

// hdf5
#define H5_USE_16_API
#include <hdf5.h>

namespace Lucee
{
/**
 * Traits class for HDF5
 */
  template <typename T> struct Hdf5Traits;

/**
 * HDF5 traits class for char
 */
  template <> struct Hdf5Traits<char> 
  {
      static hid_t hdf5Type() 
      {
        return H5T_NATIVE_CHAR;
      }
  };

/**
 * HDF5 traits class for unsigned char
 */
  template <> struct Hdf5Traits<unsigned char> 
  {
      static hid_t hdf5Type() 
      {
        return H5T_NATIVE_UCHAR;
      }
  };
  
/**
 * HDF5 traits class for short
 */
  template <> struct Hdf5Traits<short> 
  {
      static hid_t hdf5Type() 
      {
        return H5T_NATIVE_SHORT;
      }
  };

/**
 * HDF5 traits class for unsigned short
 */
  template <> struct Hdf5Traits<unsigned short> 
  {
      static hid_t hdf5Type() 
      {
        return H5T_NATIVE_USHORT;
      }
  };

  template <> struct Hdf5Traits<int> 
  {
      static hid_t hdf5Type() 
      {
        return H5T_NATIVE_INT;
      }
  };

  template <> struct Hdf5Traits<unsigned int> 
  {
      static hid_t hdf5Type() 
      {
        return H5T_NATIVE_UINT;
      }
  };

  template <> struct Hdf5Traits<long> 
  {
      static hid_t hdf5Type() 
      {
        return H5T_NATIVE_LONG;
      }
  };

  template <> struct Hdf5Traits<unsigned long> 
  {
      static hid_t hdf5Type() 
      {
        return H5T_NATIVE_ULONG;
      }
  };

  template <> struct Hdf5Traits<float> 
  {
      static hid_t hdf5Type() 
      {
        return H5T_NATIVE_FLOAT;
      }
  };

  template <> struct Hdf5Traits<double> 
  {
      static hid_t hdf5Type() 
      {
        return H5T_NATIVE_DOUBLE;
      }
  };

  template <> struct Hdf5Traits<long double>  
  {
      static hid_t hdf5Type() 
      {
        return H5T_NATIVE_LDOUBLE;
      }
  };

  template <> struct Hdf5Traits<long long> 
  {
      static hid_t hdf5Type() 
      {
        return H5T_NATIVE_LLONG;
      }
  };
}

#endif // LC_HDF5_TRAITS_H
