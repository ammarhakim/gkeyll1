/**
 * @file	lcparselfmsg.h
 *
 * @brief	Class for messaging in serial.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_PAR_SELF_MSG_H
#define LC_PAR_SELF_MSG_H

// lib includes
#include <lcdatatypes.h>
#include <lcparmsgbase.h>
#include <lcparselftmpl.h>

namespace Lucee
{
  class ParSelfMsg : public Lucee::ParMsgBase
  {
    public:
/**
 * Construct a new messenger which works in serial.
 */
      ParSelfMsg();

/**
 * Rank of process
 * 
 * @return this rank
 */
      int rank() const {
        return 0;
      }

/**
 * Number of processes taking part in messaging
 *
 * @return num of processes
 */
      unsigned numProcs() const {
        return 1;
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
  };
}

#endif // LC_PAR_SELF_MSG_H
