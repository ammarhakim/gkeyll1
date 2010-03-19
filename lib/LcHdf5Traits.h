/**
 * @file	LcHdf5Traits.h
 *
 * @brief	Class to define traits for types supported by HDF5.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
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

/** HDF5 traits class */
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

/** HDF5 traits class */
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
  
/** HDF5 traits class */
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

/** HDF5 traits class */
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

/** HDF5 traits class */
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

/** HDF5 traits class */
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

/** HDF5 traits class */
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

/** HDF5 traits class */
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

/** HDF5 traits class */
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

/** HDF5 traits class */
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

/** HDF5 traits class */
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

/** HDF5 traits class */
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
