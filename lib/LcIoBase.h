/**
 * @file	LcIoBase.h
 *
 * @brief	Base class for I/O using hierachical file system.
 */

#ifndef LC_IO_BASE_H
#define LC_IO_BASE_H

// lib includes
#include <LcDataTypes.h>
#include <LcExcept.h>
#include <LcIoNodeType.h>
#include <LcIoTmpl.h>

// etc includes
#include <loki/HierarchyGenerators.h>

// std includes
#include <vector>

namespace Lucee
{
/**
 * Provides an abstract interface for access to hierachical datasets
 */
  class IoBase 
  {
    public:
/**
 * Destory I/O object.
 */
      virtual ~IoBase();

/**
 * Create a new file for I/O.
 *
 * @param fileName the name for the file.  Assumed to be rw.
 * @return node for the file.
 */
      virtual Lucee::IoNodeType createFile(const std::string& fileName) = 0;

/**
 * Open an existing file for I/O.
 *
 * @param fileName the name for the file.
 * @param perms: the read and write permissions.  "r" or "rw".
 * @return node for the file.
 */
      virtual Lucee::IoNodeType openFile(const std::string& fileName, 
        const std::string& perms) = 0;

/**
 * Close a file node.
 *
 * @param fileNode node of file to close.
 */
      virtual void closeFile(Lucee::IoNodeType fileNode) = 0;

/**
 * Create an empty group.
 *
 * @param node the node to write under.
 * @param dataName the name of the data.
 *
 * @return the node to the written data.
 */
      virtual Lucee::IoNodeType createGroup(Lucee::IoNodeType node,
        const std::string& dataName) const = 0;

/**
 * Open a group.
 *
 * @param node the node to look in.
 * @param dataName the name of the data to open.
 *
 * @return the node to the read data.
 */
      virtual Lucee::IoNodeType openGroup(Lucee::IoNodeType node,
        const std::string& dataName) const = 0;

/**
 * Create an empty node.
 *
 * @param node the node to write under.
 * @param dataName the name of the data.
 *
 * @return the node to the written data.  This is not closed.
 */
      virtual Lucee::IoNodeType createDataSet(Lucee::IoNodeType node,
        const std::string& dataName) const = 0;

/**
 * Close a data set.
 *
 * @param node the node to be closed
 */
      virtual void closeDataSet(Lucee::IoNodeType node) const = 0;

/**
 * Open a node
 *
 * @param node the node to look in
 * @param dataName the name of the node to open
 *
 * @return the node to the read data
 */
      virtual Lucee::IoNodeType openDataSet(Lucee::IoNodeType node, 
        const std::string& dataName) const = 0;

/**
 * Write a new data set under a node.
 *
 * @param node the node to write under.
 * @param dataName the name of the data.
 * @param dataSetSize vector of the sizes of the entire data set.
 * @param dataSetBeg first index for each direction.
 * @param dataSetLen length of data to be written for each direction.
 * @param data the data to be written.
 * @return the node to the written data. The node is not closed.
 */
      template <class DATATYPE> 
      Lucee::IoNodeType
      writeDataSet(Lucee::IoNodeType node, const std::string& dataName, 
        const std::vector<size_t>& dataSetSize, const std::vector<size_t>& dataSetBeg, 
        const std::vector<size_t>& dataSetLen, const DATATYPE* data) 
      {
        return getIoPtr<DATATYPE>()->writeDataSet(node, dataName, dataSetSize,
          dataSetBeg, dataSetLen, data);
      }

/**
 * Read a data set under a node.
 *
 * @param node the node to write under.
 * @param dataName the name of the data.
 * @param dataSetBeg first index for each direction.
 * @param dataSetLen length of data to be written for each direction.
 * @param data the data to be read.
 * @return the node to the read data. The node is not closed.
 */
      template <class DATATYPE> 
      Lucee::IoNodeType
      readDataSet(Lucee::IoNodeType node, const std::string& dataName, 
        const std::vector<size_t>& dataSetBeg,
        const std::vector<size_t>& dataSetLen, DATATYPE* data)
      {
        return getIoPtr<DATATYPE>()->readDataSet(node, dataName, dataSetBeg,
          dataSetLen, data);
      }

/**
 * Write an attribute.
 *
 * @param node the node to write under
 * @param attribName the name of the attribute
 * @param attrib the attribute to be written
 */
      template <class DATATYPE> 
      void
      writeAttribute(Lucee::IoNodeType node,
	const std::string& attribName, const DATATYPE& attrib) 
      {
        getIoPtr<DATATYPE>()->writeAttribute(node, attribName, attrib);
      }

/**
 * Write a vector of attributes.
 *
 * @param node the node to write under
 * @param attribName the name of the attribute
 * @param attrib the vector attribute to be written
 */
      template <class DATATYPE>
      void
      writeVecAttribute(Lucee::IoNodeType node,
	const std::string& attribName, const std::vector<DATATYPE>& attrib) 
      {
        getIoPtr<DATATYPE>()->writeVecAttribute(node, attribName, attrib);
      }

/**
 * Write a string attribute
 *
 * @param node the node to write under
 * @param attribName the name of the attribute
 * @param attrib the string value of the attribute
 */
      virtual
      void
      writeStrAttribute(Lucee::IoNodeType node, 
        const std::string& attribName, const std::string& attrib) const = 0;

/**
 * Read an attribute.
 *
 * @param node the node to write under
 * @param attribName the name of the attribute
 * @param attrib the attribute to be read
 */
      template <class DATATYPE> void readAttribute(Lucee::IoNodeType node,
	const std::string& attribName, DATATYPE& attrib) 
      {
        getIoPtr<DATATYPE>()->readAttribute(node, attribName, attrib);
      }

/**
 * Read an attribute.
 *
 * @param node the node to write under
 * @param attribName the name of the attribute
 * @param attrib the vector attribute to be read
 */
      template <class DATATYPE> void readVecAttribute(Lucee::IoNodeType node,
	const std::string& attribName, std::vector<DATATYPE>& attrib) 
      {
        getIoPtr<DATATYPE>()->readVecAttribute(node, attribName, attrib);
      }

    protected:
/**
 * Default constructor - no dump, basename or suffix.
 *
 */
      IoBase();

/**
 * Constructor is protected, as this class cannot be made standalone.
 * Construct with suffix only.
 *
 * @param sfx the suffix for file names
 */
      IoBase(const std::string& sfx);

/**
 * Add a file to the list of open files.
 *
 * @param node the file node
 */
      virtual void addOpenFile(Lucee::IoNodeType node);

/**
 * Remove a file from the list of open files.
 *
 * @param node the file node
 */
      virtual void removeOpenFile(Lucee::IoNodeType node);

/**
 * Close any files currently open.
 */
      virtual void closeOpenFiles();

/**
 * Add a new io object. The derived class should call this to setup
 * IoBase properly.
 */
      template <typename T>
      void addIo(const IoTmpl<T>* b) 
      {
        Loki::Field<T>(ioTypeMap).ioPtr = b;
      }

    private:

/** Private copy constructor to prevent use */
      IoBase(const IoBase&);

/** Private assignment to prevent use */
      IoBase& operator=(const IoBase&);

/**
 * Get an io object of the needed type
 */
      template <typename T>
      const IoTmpl<T>* getIoPtr() {
        const IoTmpl<T>* r = Loki::Field<T>(ioTypeMap).ioPtr;
        if (r) return r;
        Lucee::Except ex;
        ex << "I/O type not set properly";
        throw ex;
      }

/** List of all open files */
      std::vector<Lucee::IoNodeType> openFiles; // List of all open files

    public:
/** Container class for all io classes */
      template <typename T>
      struct IoContainer 
      {
/** Create a new I/O container object */
          IoContainer() 
            : ioPtr(0)
          {
          }

/** Destory the container */
          virtual ~IoContainer() 
          {
            delete ioPtr;
          }
/** Pointer to a derived class of IoTmpl<T> */
          const IoTmpl<T>* ioPtr;
      };

/** Type definition for map of types to I/O containers */
      typedef Loki::GenScatterHierarchy<DataTypes_t, IoContainer> IoTypeMap;
/** Map for I/O containers */
      IoTypeMap ioTypeMap;
  };
}

#endif // LC_IO_BASE_H
