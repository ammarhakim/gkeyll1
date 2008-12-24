/**
 * @file	lcparmpimsg.h
 *
 * @brief	Class for parallel messaging using MPI library.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_PAR_MPI_MSG_H
#define LC_PAR_MPI_MSG_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#ifdef HAVE_MPI
# include <mpi.h>
#endif

// lib includes
#include <lcdatatypes.h>
#include <lcparmsgbase.h>
#include <lcparmpitmpl.h>

// std includes

namespace Lucee
{
  class ParMpiMsg : public Lucee::ParMsgBase
  {
    public:
/**
 * Construct a new MPI messenger give a set of communicating
 * processors.
 *
 * @param commProcs set of communicating processors
 */
      ParMpiMsg();

/**
 * Rank of process
 * 
 * @return this rank
 */
      int rank() const {
        int r;
        MPI_Comm_rank(_comm, &r);
        return r;
      }

/**
 * Number of processes taking part in messaging
 *
 * @return num of processes
 */
      unsigned numProcs() const {
        int np;
        MPI_Comm_size(_comm, &np);
        return np;
      }

/**
 * Split communicator into a child communicator
 *
 * @param ranks list of processors in old communicator which are to be
 * in the new communicator
 * @return new communicator
 */
      ParMsgBase* createSubComm(const std::vector<int>& ranks);


/**
 * Block till all processes hit this barrier
 */
      void barrier() const {
        MPI_Barrier(_comm);
      }

/**
 * Finish the receive started by startRecv and return a pointer to the
 * data recieved. The calling function owns the pointer and is
 * resposible for freeing it.
 *
 * @param ms message status object returned by startRecv
 */
      void * finishRecv(MsgStatus ms);

/**
 * Check status of recieve started by a startRecv. If this call
 * returns true then finishRecv can be called to immediately recieve
 * the data.
 *
 * @param ms message status object returned by startRecv
 * @return true it recieve has been completed, false otherwise
 */
      bool checkRecv(MsgStatus ms);

    private:
/**
 * Private ctor: for internal use only
 */
      ParMpiMsg(ParMpiMsg *parent, MPI_Comm comm);

      MPI_Comm _comm;
  };
}

#endif // LC_PAR_MPI_MSG_H
