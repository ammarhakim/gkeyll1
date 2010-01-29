/**
 * @file	LcIoNode.cpp
 *
 * @brief	Class representing I/O node for HDF5 file
 *
 * @version	$Id: LcIoNode.cpp 128 2009-08-20 17:15:56Z a.hakim777 $
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

// lucee includes
#include <LcIoNode.h>

namespace Lucee
{
  IoNode_v::IoNode_v(long node, IoNode_v* parent)
    : node(node), parent(parent), nodeType(1) 
  {
    if (parent)
      parent->children.push_back(this);
  }

  IoNode_v::~IoNode_v() 
  {
// loop over children killing them
    unsigned n = children.size();
    for (unsigned i=0; i<n; ++i) 
    {
      IoNode_v *ptr = children.back();
      delete ptr;
    }
    if (parent) 
    {
      if (nodeType == 1) 
      {
        H5Dclose(node);
      } 
      else if (nodeType == 0) 
      {
        H5Gclose(node);
      }
// remove us from parent
      std::vector<IoNode_v*>::iterator itr;
      for (itr = parent->children.begin(); itr != parent->children.end(); ++itr) 
      {
        if (**itr == *this) 
        {
          parent->children.erase(itr);
          break;
        }
      }
    }
    else 
    {
// we are file node
      H5Fclose(node);
    }
  }

  void 
  IoNode_v::setNodeType(unsigned type) 
  {
    nodeType = type;
  }

  bool 
  IoNode_v::operator==(const IoNode_v& v) 
  {
    return node == static_cast<const IoNode_v&>(v).node;
  }
}
