/**
 * @file	LcHdf5FileIo.h
 *
 * @brief	Class to perform I/O from HDF5 files
 *
 * @version	$Id: LcHdf5FileIo.h 128 2009-08-20 17:15:56Z a.hakim777 $
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_HDF5_FILE_IO_H
#define LC_HDF5_FILE_IO_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcExcept.h>
#include <LcIoNode.h>

// Loki includes
#include <loki/HierarchyGenerators.h>

#ifdef HAVE_MPI
#include <mpi.h>
#else
/** Set to void* */
#define MPI_Comm void*
/** Set to void* */
#define MPI_Info void*
#endif

// std includes
#include <string>

namespace Lucee
{
/**
 * HDF5 file I/O class. Provides methods for creating/opening files,
 * groups, datasets and attributes. Works for writing arbitrary types.
 */
  class Hdf5FileIo
  {
    public:
/**
 * Create a new HDF5FileIo object to perform I/O.
 *
 * @param mc the MPI communicator
 * @param mi info for the MPI communicator
 */
      Hdf5FileIo(MPI_Comm mc, MPI_Info mi);

/**
 * Create a file.
 *
 * @param fileName the name for the file.  Assumed to be rw.
 * @return node for the file
 */
      IoNode createFile(const std::string& fileName);

/**
 * Open a file.
 *
 * @param fileName the name for the file
 * @param perms: the read and write permissions.  "r" or "rw"
 * @return node for the file
 */
      IoNode openFile(const std::string& fileName, const std::string& perms);

/**
 * Close a file node
 */
      void closeFile(IoNode fileNode);

/**
 * Create an empty group
 *
 * @param node Node to create a group under
 * @param grp Name of the group
 * @return group node.
 */
      IoNode createGroup(IoNode node, const std::string& grp) const;

/**
 * Open a group
 *
 * @param node Node to look in
 * @param grp Name of the group to open
 * @return group node
 */
      IoNode openGroup(IoNode node, const std::string& dataName) const;

/**
 * Create an empty data node.
 *
 * @param node the node to write under
 * @param dataName the name of the data
 *
 * @return the node to the written data.  This is not closed.
 */
      IoNode createDataSet(IoNode node, const std::string& dataName) const;

/**
 * Open a node
 *
 * @param node the node to look in
 * @param dataName the name of the node to open
 *
 * @return the node to the read data
 */
      IoNode openDataSet(IoNode node, const std::string& dataName) const;

    private:
/** MPI communicator to use */
    MPI_Comm mpiComm;
/** MPI info object to use */
    MPI_Info mpiInfo;
  };
}

#endif // LC_HDF5_FILE_IO_H
