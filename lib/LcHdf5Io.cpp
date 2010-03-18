/**
 * @file	LcHdf5Io.cc
 *
 * @brief	Class for I/O using HDF5.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// lucee includes
#include <LcHdf5Io.h>
#include <LcHdf5IoTmpl.h>
#include <LcExcept.h>

#ifndef H5_USE_16_API
# define H5_USE_16_API
#endif
#include <hdf5.h>

namespace Lucee
{

  Hdf5Io::Hdf5Io(MPI_Comm mc, MPI_Info mi) 
  {
// Store mpi stuff, create templated IO objects
    setup(mc, mi);
  }

  Hdf5Io::~Hdf5Io() 
  {
// must close any open files.
    closeOpenFiles();
  }

  void Hdf5Io::setup(MPI_Comm mc, MPI_Info mi) 
  {
// Store mpi stuff
    mpiComm = mc;
    mpiInfo = mi;

// Create the templated io objects
    this->addIo( new Hdf5IoTmpl<char>());
    this->addIo( new Hdf5IoTmpl<unsigned char>() );
    this->addIo( new Hdf5IoTmpl<short>() );
    this->addIo( new Hdf5IoTmpl<unsigned short>() );
    this->addIo( new Hdf5IoTmpl<int>() );
    this->addIo( new Hdf5IoTmpl<unsigned int>() );
    this->addIo( new Hdf5IoTmpl<long>() );
    this->addIo( new Hdf5IoTmpl<unsigned long>() );
    this->addIo( new Hdf5IoTmpl<long long>() );
    this->addIo( new Hdf5IoTmpl<float>() );
    this->addIo( new Hdf5IoTmpl<double>() );
    this->addIo( new Hdf5IoTmpl<long double>() );
  }

  Lucee::IoNodeType
  Hdf5Io::createFile(const std::string& fileName) 
  {
// Determine access properties
    hid_t plistId = H5Pcreate(H5P_FILE_ACCESS);
#ifdef HAVE_MPI
    H5Pset_fapl_mpio(plistId, mpiComm, mpiInfo);
#endif
    long fn = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC,
      H5P_DEFAULT, plistId);
    Lucee::IoNodeTypev *fileNode = new Hdf5NodeTypev(fn);
    H5Pclose(plistId);
    if (fn < 0) 
    {
      Lucee::Except ex("Hdf5Io::Hdf5Io: unable to create file ");
      ex << "'" << fileName << "'";
      throw ex;
    }
    addOpenFile(fileNode);
    return fileNode;
  }

  Lucee::IoNodeType
  Hdf5Io::openFile(const std::string& fileName, const std::string& perms)
  {
// Determine MPI access properties
    hid_t plistId = H5Pcreate(H5P_FILE_ACCESS);
#ifdef HAVE_MPI
    H5Pset_fapl_mpio(plistId, mpiComm, mpiInfo);
#endif
// Determine whether writable
    long fn;

    if (perms.find_first_of('w') != std::string::npos) 
    {
      fn = H5Fopen(fileName.c_str(), H5F_ACC_TRUNC, plistId);
    }
    else 
    {
      fn = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, plistId);
    }
    Lucee::IoNodeTypev *fileNode = new Hdf5NodeTypev(fn);
    H5Pclose(plistId);
    if (fileNode < 0) 
    {
      Lucee::Except ex("Hdf5Io::Hdf5Io: unable to create file ");
      ex << "'" << fileName << "'";
      throw ex;
    }
    addOpenFile(fileNode);
    return fileNode;
  }

  Lucee::IoNodeType
  Hdf5Io::createGroup(Lucee::IoNodeType node, const std::string& dataName) const 
  {
    long dn = H5Gcreate(
        static_cast<Hdf5NodeTypev*>(node)->node,
        dataName.c_str(),
        0);
    Lucee::IoNodeTypev *dataNode = new Hdf5NodeTypev(dn, static_cast<Hdf5NodeTypev*>(node));
    static_cast<Hdf5NodeTypev*>(dataNode)->setNodeType(0);
    return dataNode;
  }

  Lucee::IoNodeType
  Hdf5Io::openGroup(Lucee::IoNodeType node, const std::string& dataName) const 
  {
    long dn = H5Gopen(
        static_cast<Hdf5NodeTypev*>(node)->node,
        dataName.c_str());
    Lucee::IoNodeTypev* dataNode = new Hdf5NodeTypev(dn, static_cast<Hdf5NodeTypev*>(node));
    static_cast<Hdf5NodeTypev*>(dataNode)->setNodeType(0);

    return dataNode;
  }

  Lucee::IoNodeType
  Hdf5Io::createDataSet(Lucee::IoNodeType node, const std::string& dataName) const 
  {
    hid_t fileSpace = H5Screate(H5S_SCALAR);
    long dn = H5Dcreate(
        static_cast<Hdf5NodeTypev*>(node)->node,
        dataName.c_str(),
        H5T_NATIVE_INT,
        fileSpace,
        H5P_DEFAULT);
    Lucee::IoNodeTypev *dataNode = new Hdf5NodeTypev(dn, static_cast<Hdf5NodeTypev*>(node));
    H5Sclose(fileSpace);

    return dataNode;
  }

  Lucee::IoNodeType
  Hdf5Io::openDataSet(Lucee::IoNodeType node, const std::string& dataName) const 
  {
    long dn = H5Dopen(
        static_cast<Hdf5NodeTypev*>(node)->node,
        dataName.c_str());
    Lucee::IoNodeTypev* dataNode = new Hdf5NodeTypev(dn, static_cast<Hdf5NodeTypev*>(node));

    return dataNode;
  }

  void Hdf5Io::closeFile(Lucee::IoNodeTypev *fileNode) 
  {
    removeOpenFile(fileNode); // Remove reference in base class
    delete fileNode;
  }

  void Hdf5Io::closeDataSet(Lucee::IoNodeTypev *node) const 
  {
    delete node;
  }

  void
  Hdf5Io::writeStrAttribute(Lucee::IoNodeType node, 
    const std::string& attribName, const std::string& attrib) const 
  {
// Add attribute that saves a vector of chars
    hsize_t size = attrib.size();
    hid_t attrDS = H5Screate(H5S_SCALAR);
// Create attribute to store the number of elements
    hid_t attrTP = H5Tcopy(H5T_C_S1);
    H5Tset_size(attrTP, size + 1);
    hid_t attrID = H5Acreate(static_cast<Hdf5NodeTypev*>(node)->node,
      attribName.c_str(), attrTP, attrDS, H5P_DEFAULT);
// Write out the attribute
    H5Awrite(attrID, attrTP, attrib.c_str());
//this->comm()->barrier();
// Close the attribute
    H5Aclose(attrID);
// Close the attribute type
    H5Tclose(attrTP);
// Close the attribute dataspace
    H5Sclose(attrDS);
  }
}
