/**
 * @file	lcparmsgtmpl.h
 *
 * @brief	Base class for parallel message passing classes.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_PAR_MSG_OPS_H
#define LC_PAR_MSG_OPS_H

namespace Lucee
{

// list of supported operations for all-reduce operation
  enum ParMsgOps
  {
    PAR_MSG_NOP,
    PAR_MSG_MIN,
    PAR_MSG_MAX,
    PAR_MSG_AND,
  };

/**
 * Generic reduction operation wrapper
 */
  class ParMsgOp
  {
    protected:
/**
 * Constructor
 */
      ParMsgOp() 
      {
      }
  };
}

#endif // LC_PAR_MSG_OPS_H
