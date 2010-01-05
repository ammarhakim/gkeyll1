/**
 * @file	LcHdf5File.h
 *
 * @brief	Class to represent HDF5 file.
 *
 * @version	$Id: LcHdf5File.h 128 2009-08-20 17:15:56Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_HDF5_FILE_H
#define LC_HDF5_FILE_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#else
/** Set to void* */
#define MPI_Comm void*
/** Set to void* */
#define MPI_Info void*
#endif

// HDF5 includes
#include <hdf5.h>

// std includes
#include <string>

namespace Lucee
{
/**
 * 
 */
  class Hdf5File
  {
    public:
/**
 * Open a HDf5 file for I/O. The file can be opened either for reading
 * or reading/writing.
 *
 * @param fname Name of file to open.
 * @param perms Either "r" or "rw".
 */
      Hdf5File(const std::string& fname, const std::string& perms);

    private:
/** Name of HDF5 file */
      std::string fileName;
  };
}

#endif //  LC_HDF5_FILE_H
