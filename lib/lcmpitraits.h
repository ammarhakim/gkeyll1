/**
 * @file	lcmpitraits.h
 *
 * @brief	Class to define traits for types supported by MPI.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_MPI_TRAITS_H
#define LC_MPI_TRAITS_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#ifdef HAVE_MPI
# include <mpi.h>
#endif

#include <lcparmsgbase.h>

namespace Lucee
{
/**
 * Type traits for use in MPI messengers
 */
  template<typename T> class MpiTraits;

// char
  template<> class MpiTraits<char> 
  {
    public:
      static MPI_Datatype mpiType() 
      {
        return MPI_CHAR;
      };
  };
// unsigned char
  template<> class MpiTraits<unsigned char> 
  { 
    public:
      static MPI_Datatype mpiType() 
      {
        return MPI_UNSIGNED_CHAR;
      }; 
  };
// short
  template<> class MpiTraits<short> 
  {
    public:
      static MPI_Datatype mpiType() 
      {
        return MPI_SHORT;
      }; 
  };
// unsigned short
  template<> class MpiTraits<unsigned short> 
  {
    public:
      static MPI_Datatype mpiType() 
      {
        return MPI_UNSIGNED_SHORT;
      }; 
  };
// int
  template<> class MpiTraits<int> 
  {
    public:
      static MPI_Datatype mpiType() 
      {
        return MPI_INT;
      };
  };
// unsigned
  template<> class MpiTraits<unsigned> 
  {
    public:
      static MPI_Datatype mpiType() 
      {
        return MPI_UNSIGNED;
      };
  };
// long
  template<> class MpiTraits<long> 
  {
    public:
      static MPI_Datatype mpiType() 
      {
        return MPI_LONG;
      };
  };
// unsigned long
  template<> class MpiTraits<unsigned long> 
  {
    public:
      static MPI_Datatype mpiType() 
      {
        return MPI_UNSIGNED_LONG;
      };
  };
// float
  template<> class MpiTraits<float> 
  {
    public:
      static MPI_Datatype mpiType() 
      {
        return MPI_FLOAT;
      };
  };
// double
  template<> class MpiTraits<double> 
  {
    public:
      static MPI_Datatype mpiType() 
      {
        return MPI_DOUBLE;
      };
  };
// long double
  template<> class MpiTraits<long double> 
  {
    public:
      static MPI_Datatype mpiType() 
      {
        return MPI_LONG_DOUBLE;
      };
  };
// long long int
  template<> class MpiTraits<long long int> 
  {
    public:
      static MPI_Datatype mpiType() 
      {
        return MPI_LONG_LONG_INT;
      };
  };

/**
 * Type traits for use in MPI all-reduce operators
 */
  template<const unsigned T> class MpiAllReduce;

// MIN
  template<> class MpiAllReduce<PAR_MSG_MIN> 
  { 
    public:
      static MPI_Op mpiOp() 
      {
        return MPI_MIN;
      };
  };
// MAX
  template<> class MpiAllReduce<PAR_MSG_MAX> 
  { 
    public:
      static MPI_Op mpiOp() 
      {
        return MPI_MAX;
      };
  };
}

#endif // LC_MPI_TRAITS_H
