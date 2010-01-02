/**
 * @file	LcHdf5FileIo.cc
 *
 * @brief	Class to perform I/O from HDF5 files
 *
 * @version	$Id: LcHdf5FileIo.cc 128 2009-08-20 17:15:56Z a.hakim777 $
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcHdf5FileIo.h>
#include <LcExcept.h>

namespace Lucee

{
  Hdf5FileIo::Hdf5FileIo(MPI_Comm mc, MPI_Info mi)
    : mpiComm(mc), mpiInfo(mi)
  {
  }

  IoNode
  Hdf5FileIo::createFile(const std::string& fileName)
  {
// Determine access properties
    hid_t plistId = H5Pcreate(H5P_FILE_ACCESS);
#ifdef HAVE_MPI
    H5Pset_fapl_mpio(plistId, mpiComm, mpiInfo);
#endif
    long fn = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC,
      H5P_DEFAULT, plistId);
    IoNode_v *fileNode = new IoNode_v(fn);
    H5Pclose(plistId);
    if (fn < 0) {
      Except lcex("Hdf5FileIo::createFile: unable to create file ");
      lcex << "'" << fileName << "'";
      throw lcex;
    }
    addOpenFile(fileNode);
    return fileNode;
  }    

  IoNode
  Hdf5FileIo::openFile(const std::string& fileName, const std::string& perms)
  {
// Determine MPI access properties
    hid_t plistId = H5Pcreate(H5P_FILE_ACCESS);
#ifdef HAVE_MPI
    H5Pset_fapl_mpio(plistId, mpiComm, mpiInfo);
#endif
// Determine whether writable
    long fn;

    if (perms.find_first_of('w') != std::string::npos) {
      fn = H5Fopen(fileName.c_str(), H5F_ACC_TRUNC, plistId);
    }
    else {
      fn = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, plistId);
    }
    IoNode_v *fileNode = new IoNode_v(fn);
    H5Pclose(plistId);
    if (fileNode < 0) {
      Except lcex("Hdf5FileIo:::openFile unable to open file ");
      lcex << "'" << fileName << "'";
      throw lcex;
    }
    addOpenFile(fileNode);
    return fileNode;
  }

  void
  Hdf5FileIo::closeFile(IoNode fileNode)
  {
  }

  IoNode
  Hdf5FileIo::createGroup(IoNode node, const std::string& grp) const
  {
    return node;
  }
  
  IoNode
  Hdf5FileIo::openGroup(IoNode node, const std::string& dataName) const
  {
    return node;
  }

  IoNode
  Hdf5FileIo::createDataSet(IoNode node, const std::string& dataName) const
  {
    return node;
  }
  
  IoNode
  Hdf5FileIo::openDataSet(IoNode node, const std::string& dataName) const
  {
    return node;
  }

}
