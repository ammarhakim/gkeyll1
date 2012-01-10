/**
 * @file	LcNamedDataNode.h
 *
 * @brief	Node to hold named data.
 */

#ifndef LC_NAMED_DATA_NODE_H
#define LC_NAMED_DATA_NODE_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcArray.h>
#include <LcVector.h>

namespace Lucee
{
/**
 * This class holds data of different types and forms a node in a
 * tree. The data can be accessed using names. Nodes are usually
 * initialized from data stored in files.
 */
  class NamedDataNode
  {
    public:

    private:
/** Private structure to hold map of name -> data */
      template <typename T>
      struct DataMap
      {
/** Map storing name -> data */
          std::map<std::string, T> dataMap;
      };
  };
}

#endif // LC_NAMED_DATA_NODE_H
