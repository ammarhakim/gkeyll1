/**
 * @file	lciotmpl.h
 *
 * @brief	Base class for I/O nodes.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_IO_TMPL_H
#define LC_IO_TMPL_H

// lib includes
#include <lcionodetype.h>

// std includes
#include <string>
#include <vector>

namespace Lucee
{
/**
 * IoTmpl is the base class for access to a hierarchical file
 * system with groups, data sets, and attributes for those datasets.
 * The exemplar is HDF5, but one may eventually other systems, like
 * netCDF, PDB, or ?
 *
 * This is the base class.  It will provide for output of the
 * basic data structures of the STL library and boost objects.
 * Derived classes will allow for output of Facets objects, such
 * as Facets arrays, whether distributed or not.
 */
  template <class DATATYPE>
  class IoTmpl 
  {
    public:
/**
 * Virtual destructor
 */
      virtual ~IoTmpl() 
      {
      }

/**
 * Write a new data set under a node.
 *
 * @param node the node to write under.
 * @param dataName the name of the data.
 * @param dataSetSize vector of the sizes of the entire data set.
 * @param dataSetBeg first index for each direction.
 * @param dataSetLen length of data to be written for each direction.
 * @param data the data to be written.
 *
 * @return the node to the written data.
 */
      virtual Lucee::IoNodeType writeDataSet(Lucee::IoNodeType node,
        const std::string& dataName, const std::vector<size_t>& dataSetSize,
        const std::vector<size_t>& dataSetBeg,
        const std::vector<size_t>& dataSetLen, const DATATYPE* data) const = 0;

/**
 * Read a new data set under a node
 *
 * @param node the node to write under
 * @param dataName the name of the data
 * @param dataSetBeg first index for each direction
 * @param dataSetLen length of data to be written for each direction
 * @param data the data to be written
 *
 * @return the node to the written data
 */
      virtual Lucee::IoNodeType readDataSet(Lucee::IoNodeType node,
        const std::string& dataName, const std::vector<size_t>& dataSetBeg,
        const std::vector<size_t>& dataSetLen, DATATYPE* data) const = 0;

/**
 * Write an attribute.
 *
 * @param node the node to write under
 * @param attribName the name of the attribute
 * @param attrib the attribute to be written
 */
      virtual void writeAttribute(Lucee::IoNodeType node,
        const std::string& attribName, const DATATYPE& attrib) const = 0;

/**
 * Write an attribute.
 *
 * @param node the node to write under
 * @param attribName the name of the attribute
 * @param attrib the vector attribute to be written
 */
      virtual void writeVecAttribute(Lucee::IoNodeType node,
        const std::string& attribName, const std::vector<DATATYPE>& attrib) const = 0;

/**
 * Read an attribute.
 *
 * @param node the node to which this attribute belongs
 * @param attribName the name of the attribute
 * @param attrib the attribute to be read
 */
      virtual void readAttribute(Lucee::IoNodeType node,
        const std::string& attribName, DATATYPE& attrib) const = 0;

/**
 * Read an attribute.
 *
 * @param node the node to which this attribute belongs
 * @param attribName the name of the attribute
 * @param attrib the vector attribute to be read
 */
      virtual void readVecAttribute(Lucee::IoNodeType node,
        const std::string& attribName, std::vector<DATATYPE>& attrib) const = 0;

    protected:
/**
 * Constructor is protected, as this class cannot be made standalone.
 */
      IoTmpl() 
      {
      }

    private:
/** Private copy constructor to prevent use */
      IoTmpl(const IoTmpl<DATATYPE>&);

/** Private assignment to prevent use */
      IoTmpl<DATATYPE>& operator=(const IoTmpl<DATATYPE>&);
  };

}

#endif // LC_IO_TMPL_H

