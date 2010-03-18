/**
 * @file	LcIoNodeType.h
 *
 * @brief	Base class for I/O nodes.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_IO_NODE_TYPE_H
#define LC_IO_NODE_TYPE_H

namespace Lucee
{

/**
 * Provides a means for derived messengers to return implimentation
 * specific message status flags and data. The returned type is an
 * opaque pointer and should not be directly fiddled around with.
 */
  struct IoNodeTypev 
  {
/** Destructor */
      virtual ~IoNodeTypev() 
      {
      }

/**
 * A comparison operator.
 *
 * @param v the node to which we are comparing ourself.
 * @return a boolean specifying whether or not we are equal.
 */
      virtual bool operator==(const IoNodeTypev& v) = 0;
  };

/**
 * A simple typedef to easily refer to an IoNodeType pointer w/o
 * having to know about it
 */
  typedef IoNodeTypev* IoNodeType;
}

#endif // LC_IO_NODE_TYPE_H
