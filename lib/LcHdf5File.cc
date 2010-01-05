/**
 * @file	LcHdf5File.cc
 *
 * @brief	Class to represent HDF5 file.
 *
 * @version	$Id: LcHdf5File.cc 128 2009-08-20 17:15:56Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

// lLucee include
#include <LcHdf5File.h>

namespace Lucee
{
  Hdf5File::Hdf5File(const std::string& fname, const std::string& perms)
    : fileName(fname)
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
}
