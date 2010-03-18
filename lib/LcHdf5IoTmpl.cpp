/**
 * @file	LcHdf5IoTmpl.cpp
 *
 * @brief	HDF5 I/O class.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcExcept.h>
#include <LcHdf5IoTmpl.h>
#include <LcHdf5Traits.h>

// std includes
#include <iostream>
#include <string>

/**  Use HDF5 1.6 API */
#define H5_USE_16_API
#include <hdf5.h>
/**  Use HDF5 1.6 API */
#define HDF5_FAIL -1

namespace Lucee
{
  template <class DATATYPE>
  Hdf5IoTmpl<DATATYPE>::Hdf5IoTmpl() 
  {
  }

  template <class DATATYPE>
  Hdf5IoTmpl<DATATYPE>::~Hdf5IoTmpl() 
  {
  }

  template <class DATATYPE>
  Lucee::IoNodeType
  Hdf5IoTmpl<DATATYPE>::writeDataSet(Lucee::IoNodeType node,
    const std::string& dataName,
    const std::vector<size_t>& dataSetSize,
    const std::vector<size_t>& dataSetBeg,
    const std::vector<size_t>& dataSetLen,
    const DATATYPE* data) const 
  {
    // Create data shape, then data
    hsize_t* dataShape = new hsize_t[dataSetSize.size()];
    for (size_t i=0; i<dataSetSize.size(); ++i) dataShape[i] = dataSetSize[i];
    hid_t fileSpace = H5Screate_simple(dataSetSize.size(), dataShape, NULL);
    long dn = H5Dcreate(static_cast<Hdf5NodeTypev*>(node)->node,
      dataName.c_str(),
      Hdf5Traits<DATATYPE>::hdf5Type(),
      fileSpace,
      H5P_DEFAULT);
    Lucee::IoNodeTypev *dataNode =
      new Hdf5NodeTypev(dn, static_cast<Hdf5NodeTypev*>(node));

    H5Sclose(fileSpace);
    delete [] dataShape;
    if (dn == HDF5_FAIL) 
    {
      throw Lucee::Except("Hdf5IoTmpl::writeDataSet: H5Dcreate failed.");
    }

    // Create the data region to be written
#ifdef NEW_H5S_SELECT_HYPERSLAB_IFC
    hsize_t* dataStart = new hsize_t[dataSetSize.size()];
#else
    hssize_t* dataStart = new hssize_t[dataSetSize.size()];
#endif
    for (size_t i=0; i<dataSetBeg.size(); ++i) dataStart[i] = dataSetBeg[i];
    for (size_t i=dataSetBeg.size(); i<dataSetSize.size(); ++i) dataStart[i] = 0;
    hsize_t* dataExtent = new hsize_t[dataSetSize.size()];
    for (size_t i=0; i<dataSetLen.size(); ++i) dataExtent[i] = dataSetLen[i];
    for (size_t i=dataSetLen.size(); i<dataSetSize.size(); ++i) dataExtent[i] = 1;
    hid_t memSpace = H5Screate_simple(dataSetSize.size(), dataExtent, NULL);
    hid_t dataSpace = H5Dget_space(dn);
    herr_t ret = H5Sselect_hyperslab(dataSpace, H5S_SELECT_SET, dataStart,
      NULL, dataExtent, NULL);
    delete [] dataStart;
    delete [] dataExtent;
    if (ret == HDF5_FAIL) 
    {
      H5Sclose(memSpace);
      H5Dclose(dataSpace);
      throw Lucee::Except("Hdf5IoTmpl::writeDataSet: H5Sselect_hyperslab failed.");
    }

    // Set for collective I/O
#ifdef HAVE_MPI
    hid_t xferPropList = H5Pcreate(H5P_DATASET_XFER);
    ret = H5Pset_dxpl_mpio(xferPropList, H5FD_MPIO_COLLECTIVE);
    if (ret == HDF5_FAIL) 
    {
      H5Sclose(memSpace);
      H5Dclose(dataSpace);
      throw Lucee::Except("Hdf5IoTmpl::writeDataSet: H5Pset_dxpl_mpio failed.");
    }
#else
    hid_t xferPropList = H5P_DEFAULT;
#endif

    // Write the data
    ret = H5Dwrite(dn, Hdf5Traits<DATATYPE>::hdf5Type(),
      memSpace, dataSpace, xferPropList, data);

    // Close all HDF5 data
#ifdef HAVE_MPI
    H5Pclose(xferPropList);
#endif
    H5Sclose(memSpace);
    H5Sclose(dataSpace);

    // If failed, throw exception
    if (ret == HDF5_FAIL) 
    {
      throw Lucee::Except("Hdf5IoTmpl::writeDataSet: H5Dwrite failed.");
    }

    return dataNode;
  }

  template <class DATATYPE>
  Lucee::IoNodeType
  Hdf5IoTmpl<DATATYPE>::readDataSet(Lucee::IoNodeType node,
    const std::string& dataName,
    const std::vector<size_t>& dataSetBeg,
    const std::vector<size_t>& dataSetLen,
    DATATYPE* data) const 
  {
    // Open the data for reading
    long dn = H5Dopen(static_cast<Hdf5NodeTypev*>(node)->node,
      dataName.c_str());
    Lucee::IoNodeTypev *dataNode =
      new Hdf5NodeTypev(dn,
        static_cast<Hdf5NodeTypev*>(node));

    // Create the data space and memory space for reading
#ifdef NEW_H5S_SELECT_HYPERSLAB_IFC
    hsize_t* dataStart = new hsize_t[dataSetBeg.size()];
#else
    hssize_t* dataStart = new hssize_t[dataSetBeg.size()];
#endif
    for (size_t i=0; i<dataSetBeg.size(); ++i) dataStart[i] = dataSetBeg[i];
    hsize_t* dataExtent = new hsize_t[dataSetBeg.size()];
    for (size_t i=0; i<dataSetLen.size(); ++i) dataExtent[i] = dataSetLen[i];
    for (size_t i=dataSetLen.size(); i<dataSetBeg.size(); ++i) dataExtent[i] = 1;
    hid_t memSpace = H5Screate_simple(dataSetBeg.size(), dataExtent, NULL);
    hid_t dataSpace = H5Dget_space(dn);
    herr_t ret = H5Sselect_hyperslab(dataSpace, H5S_SELECT_SET, dataStart,
      NULL, dataExtent, NULL);
    delete [] dataStart;
    delete [] dataExtent;
    if (ret == HDF5_FAIL) 
    {
      H5Sclose(memSpace);
      H5Sclose(dataSpace);
      throw Lucee::Except("Hdf5IoTmpl::readDataSet: hyperslab selection failed.");
    }

    // Read the data
    ret = H5Dread(dn, Hdf5Traits<DATATYPE>::hdf5Type(),
      memSpace, dataSpace, H5P_DEFAULT, data);

    // Close all HDF5 data
    H5Sclose(memSpace);
    H5Sclose(dataSpace);

    // If failed, throw exception
    if (ret == HDF5_FAIL) 
    {
      throw Lucee::Except("Hdf5IoTmpl::readDataSet: H5Dread failed.");
    }

    // Return node
    return dataNode;
  }

  template <class DATATYPE>
  void
  Hdf5IoTmpl<DATATYPE>::writeAttribute(Lucee::IoNodeType node,
    const std::string& attribName,
    const DATATYPE& attrib) const 
  {
    // Add attribute that saves a scalar
    hid_t attrDS = H5Screate(H5S_SCALAR);
    // Create attribute to store the num of elements
    hid_t attrId = H5Acreate(static_cast<Hdf5NodeTypev*>(node)->node,
      attribName.c_str(),
      Hdf5Traits<DATATYPE>::hdf5Type(),
      attrDS,
      H5P_DEFAULT);
    // Write out the attribute
    H5Awrite(attrId, Hdf5Traits<DATATYPE>::hdf5Type(), &attrib);
    // Close the attribute
    H5Aclose(attrId);
    // Close the attribute dataspace
    H5Sclose(attrDS);
  }

  template <class DATATYPE>
  void
  Hdf5IoTmpl<DATATYPE>::writeVecAttribute(Lucee::IoNodeType node,
    const std::string& attribName,
    const std::vector<DATATYPE>& attrib) const 
  {
    // Add attribute that saves a vector
    hsize_t size = attrib.size();
    DATATYPE *data = new DATATYPE[attrib.size()];
    for (unsigned i = 0; i < attrib.size(); ++i) data[i] = attrib[i];
    hid_t attrDS = H5Screate_simple(1, &size, NULL);
    // Create attribute to store the num of elements
    hid_t attrId = H5Acreate(static_cast<Hdf5NodeTypev*>(node)->node,
      attribName.c_str(),
      Hdf5Traits<DATATYPE>::hdf5Type(),
      attrDS,
      H5P_DEFAULT);
    // Write out the attribute
    H5Awrite(attrId, Hdf5Traits<DATATYPE>::hdf5Type(), data);
    // Close the attribute
    H5Aclose(attrId);
    // Close the attribute dataspace
    H5Sclose(attrDS);
    // Tidy up
    delete [] data;
  }

  template <class DATATYPE>
  void
  Hdf5IoTmpl<DATATYPE>::readAttribute(Lucee::IoNodeType node,
    const std::string& attribName,
    DATATYPE& attrib) const 
  {
    // Read a scalar attribute
    hid_t attrId = H5Aopen_name(static_cast<Hdf5NodeTypev*>(node)->node,
      attribName.c_str());
    H5Aread(attrId, Hdf5Traits<DATATYPE>::hdf5Type(), &attrib);
    // Close the attribute
    H5Aclose(attrId);
  }

  template <class DATATYPE>
  void
  Hdf5IoTmpl<DATATYPE>::readVecAttribute(Lucee::IoNodeType node,
    const std::string& attribName,
    std::vector<DATATYPE>& attrib) const 
  {
    // Read a vector attribute
    hid_t attrId = H5Aopen_name(static_cast<Hdf5NodeTypev*>(node)->node,
      attribName.c_str());
    hid_t attrSp = H5Aget_space(attrId);
    hsize_t len;
    H5Sget_simple_extent_dims(attrSp, &len, NULL);
    DATATYPE* data = new DATATYPE[len];
    H5Aread(attrId, Hdf5Traits<DATATYPE>::hdf5Type(), data);
    attrib.clear();
    for (unsigned i = 0; i < len; ++i) attrib.push_back(data[i]);
    delete [] data;
    // Close the attribute
    H5Aclose(attrId);
  }

  //
  // Instantiate only needed classes
  //
  template class Hdf5IoTmpl<char>;
  template class Hdf5IoTmpl<unsigned char>;
  template class Hdf5IoTmpl<short>;
  template class Hdf5IoTmpl<unsigned short>;
  template class Hdf5IoTmpl<int>;
  template class Hdf5IoTmpl<unsigned int>;
  template class Hdf5IoTmpl<long>;
  template class Hdf5IoTmpl<unsigned long>;
  template class Hdf5IoTmpl<long long>;
  template class Hdf5IoTmpl<float>;
  template class Hdf5IoTmpl<double>;
  template class Hdf5IoTmpl<long double>;
}
