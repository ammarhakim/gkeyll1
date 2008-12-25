/**
 * @file	lcparmsgop.h
 *
 * @brief	Base class for parallel message passing classes.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_PAR_MSG_OP_H
#define LC_PAR_MSG_OP_H

namespace Lucee
{
/** List of supported operations for all-reduce operations */
  enum ParMsgOpCode
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

#endif // LC_PAR_MSG_OP_H
