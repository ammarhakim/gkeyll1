/**
 * @file	lcparselftmpl.h
 *
 * @brief	Class for messaging in serial.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_PAR_SELF_TMPL_H
#define LC_PAR_SELF_TMPL_H

// lib includes
#include <lcparmsgtmplbase.h>

// std includes
#include <iostream>

namespace Lucee
{
/**
 * Self communicator specific message status wrapper
 */
  struct SelfMsgStatus_v : public MsgStatus_v
  {
  };

  template<typename T>
  class ParSelfTmpl : public Lucee::ParMsgTmplBase<T>
  {
    public:
/** Create a new self-communication object */
      ParSelfTmpl()
      {
      }

/**
 * Send an array to another rank.
 *
 * @param arr std::vector of data being sent.
 * @param recvRank rank that will receive array.
 * @param tag Tag to attach to message.
 */
      void send(const std::vector<T>& arr, unsigned recvRank, int tag);

/**
 * Send an array to another rank.
 * 
 * @param num number of elements to send.
 * @param arr array of length 'num' of data being sent.
 * @param recvRank rank that will receive array.
 * @param tag Tag to attach to message.
 */
      void send(unsigned num, T* arr, unsigned recvRank, int tag);

/**
 * Receive an array from another rank.
 *
 * @param num number of elements to recv.
 * @param array array that is filled with received values.
 * @param sendRank rank that sent array.
 * @param tag Tag to attach to message.
 */
      void recv(unsigned num, std::vector<T>& array, unsigned sendRank, int tag);

/**
 * Receive an array from another rank.
 *
 * @param num number of elements to receive.
 * @param array array that is filled with received values.
 * @param sendRank rank that sent array.
 * @param tag Tag to attach to message.
 */
      void recv(unsigned num, T* array, unsigned sendRank, int tag);

/**
 * Receive an array from another rank. This is a non-blocking call.
 * 
 * @param num number of elements to receive.
 * @param sendRank rank that sent array.
 * @param tag Tag to attach to message.
 * @return message status
 */
      MsgStatus startRecv(unsigned num, unsigned sendRank, int tag);

/**
 * Reduce data to all ranks
 * 
 * @param num Number of elements being reduced
 * @param sendBuff buffer to reduce
 * @param recvBuff buffer to recieve reduced data
 * @param op operation to perform. These are one of those listed in MsgOp enum
 */
      void allReduce(unsigned num, T* sendBuff, T* recvBuff, Lucee::ParMsgOpCode op);
  };
}

#endif // LC_PAR_SELF_TMPL_H

