/**
 * @file	lcmpimsgop.cc
 *
 * @brief	Class for parallel messaging using MPI library.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

// lib includes
#include <lcmpimsgop.h>

namespace Lucee
{
  MpiMsgOp::MpiMsgOp()
  {
    _opMap[PAR_MSG_MIN] = MPI_MIN;
    _opMap[PAR_MSG_MAX] = MPI_MAX;
    _opMap[PAR_MSG_AND] = MPI_LAND;
  }

  MPI_Op
  MpiMsgOp::getOp(Lucee::ParMsgOps op)
  {
    return _opMap[op];
  }
}
