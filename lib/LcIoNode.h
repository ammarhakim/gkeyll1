/**
 * @file	LcIoNode.h
 *
 * @brief	Class representing I/O node for HDF5 file
 *
 * @version	$Id: LcIoNode.h 128 2009-08-20 17:15:56Z a.hakim777 $
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_IO_NODE_H
#define LC_IO_NODE_H

// HDF5 includes
#define H5_USE_16_API
#include <hdf5.h>

// std includes
#include <vector>

namespace Lucee
{
/**
 * I/O node to represent nodes in HDF5 files. This class represents
 * files, groups and datasets. This is used as an opaque struct to
 * prevent access to internals.
 */
  struct IoNode_v 
  {
/**
 * Create a new node object to represent a HDF5 node.
 *
 * @param node the node number of this node
 * @param parent the parent of this node
 */
      IoNode_v(long node=0, IoNode_v* parent=0);

/** Destructor */
      virtual ~IoNode_v();

/**
 * Sets the type of this node
 *
 * @param type the type of this node
 */
      void setNodeType(unsigned type);

/**
 * Comparison operator
 *
 * @param v the node to which we are comparing ourself
 * @return boolean specifying whether we are the same
 */
      bool operator==(const IoNode_v& v);

/** This node's node number */
      long node;
/** This node's parent node */
      IoNode_v *parent;
/** This node's list of children */
      std::vector<IoNode_v*> children;
/** This node's type */
      int nodeType;
  };

/** Typedef to make node object opaque */
  typedef IoNode_v* IoNode;
}

#endif // LC_IO_NODE_H
