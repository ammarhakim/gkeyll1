/**
 * @file	LcHdf5Io.h
 *
 * @brief	Class for I/O using HDF5.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_HDF5_IO_H
#define LC_HDF5_IO_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcIoBase.h>

// check if we are using MPI libraries.
#ifdef HAVE_MPI
# include <mpi.h>
#else
/** Communicator set to void if not using MPI */
# define MPI_Comm void*
/** Info set to void if not using MPI */
# define MPI_Info void*
#endif

namespace Lucee
{
/**
 * Hdf5Io is the interface for the HDF5 implementation of HDF5.
 */
  class Hdf5Io : public Lucee::IoBase
  {
    public:
/**
 * Constructor creates the individual templated writers
 *
 * @param mc the MPI communicator
 * @param mi info for the MPI communicator
 */
      Hdf5Io(MPI_Comm mc, MPI_Info mi);

/**
 * Constructor creates the individual templated writers
 *
 * @param bn base name for the dump
 * @param d the dump number
 * @param mc the MPI communicator
 * @param mi info for the MPI communicator
 */
      Hdf5Io(const std::string& bn, int d, MPI_Comm mc, MPI_Info mi);

/**
 * Virtual destructor
 */
      virtual ~Hdf5Io();

/**
 * Create a file.
 *
 * @param fileName the name for the file.  Assumed to be rw.
 *
 * @return node for the file
 */
      virtual Lucee::IoNodeType createFile(const std::string& fileName);

/**
 * Open a file.
 *
 * @param fileName the name for the file
 * @param perms: the read and write permissions.  "r" or "rw"
 *
 * @return node for the file
 */
      virtual Lucee::IoNodeType openFile(const std::string& fileName,
        const std::string& perms);

/**
 * Create an empty group
 *
 * @param node the node to write under
 * @param dataName the name of the data
 *
 * @return the node to the written data
 */
      virtual Lucee::IoNodeType createGroup(Lucee::IoNodeType node, 
        const std::string& dataName) const;

/**
 * Open a group
 *
 * @param node the node to look in
 * @param dataName the name of the data to open
 *
 * @return the node to the read data
 */
      virtual Lucee::IoNodeType openGroup(Lucee::IoNodeType node, 
        const std::string& dataName) const;

/**
 * Create an empty node
 *
 * @param node the node to write under
 * @param dataName the name of the data
 *
 * @return the node to the written data
 */
      virtual Lucee::IoNodeType createDataSet(Lucee::IoNodeType node, 
        const std::string& dataName) const;

/**
 * Open a node
 *
 * @param node the node to look in
 * @param dataName the name of the node to open
 *
 * @return the node to the read data
 */
      virtual Lucee::IoNodeType openDataSet(Lucee::IoNodeType node, 
        const std::string& dataName) const;

/**
 * Get the top node = file node
 */
      virtual void closeFile(Lucee::IoNodeType fileNode);

/**
 * Close a data set.
 *
 * @param node the node to be closed
 */
      virtual void closeDataSet(Lucee::IoNodeType node) const;

/**
 * Write a string attribute
 *
 * @param node the node to write under
 * @param attribName the name of the attribute
 * @param attrib the string value of the attribute
 */
      void writeStrAttribute(Lucee::IoNodeType node, 
        const std::string& attribName, const std::string& attrib) const;

    private:
/** Private copy constructor to prevent use */
      Hdf5Io(const Hdf5Io&);

/** Private assignment to prevent use */
      Hdf5Io& operator=(const Hdf5Io&);

/**
 * Set up the templated IO objects and the mpi stuff
 */
      void setup(MPI_Comm mc, MPI_Info mi);

/** The communicator */
      MPI_Comm mpiComm;
/** The mpi info */
      MPI_Info mpiInfo;
  };
}

#endif // LC_HDF5_IO_H
