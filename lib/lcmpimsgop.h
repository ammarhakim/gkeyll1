/**
 * @file	lcmpimsgop.h
 *
 * @brief	Class for parallel messaging using MPI library.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_MPI_MSG_OP_H
#define LC_MPI_MSG_OP_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lib includes
#include <lcparmsgop.h>

#ifdef HAVE_MPI
# include <mpi.h>
#endif

// std includes
#include <map>

namespace Lucee
{
  class MpiMsgOp : public Lucee::ParMsgOp
  {
    public:
      MpiMsgOp();

/**
 * Returns an MPI operation corresponding to a string key.
 *
 * @param op an operation.
 * @return an MPI operation.
 */
      MPI_Op getOp(Lucee::ParMsgOps op);

    private:
      std::map<Lucee::ParMsgOps, MPI_Op> _opMap;
  };
}

#endif // __wxmpimsgop__
