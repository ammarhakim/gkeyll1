/**
 * @file	lcparmpitmpl.h
 *
 * @brief	Class for parallel messaging using MPI library.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_PAR_MPI_TMPL_H
#define LC_PAR_MPI_TMPL_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#ifdef HAVE_MPI
# include <mpi.h>
#endif

// lib includes
#include <lcparmsgtmplbase.h>
#include <lcmpitraits.h>
#include <lcmpimsgop.h>

// std includes
#include <iostream>

namespace Lucee
{
/**
 * Mpi specific message status wrapper
 */
  struct MpiMsgStatus_v : public Lucee::MsgStatus_v
  {
/** Create a new status object */
      MpiMsgStatus_v()
        : data(0) 
      {
      }
/** Send/recv equest object */
      MPI_Request request;
/** Pointer to memory for recv */
      void *data;
  };

  template<typename T>
  class ParMpiTmpl : public Lucee::ParMsgTmplBase<T>
  {
    public:
/** 
 * Create a new messaging object using MPI communicator.
 *
 * @param comm MPI communicator.
 */
      ParMpiTmpl(MPI_Comm comm)
        : _comm(comm), _sending(false) 
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
 * @param num number of elements to send
 * @param arr array of length 'num' of data being sent
 * @param recvRank rank that will receive array
 * @param tag Tag to attach to message.
 */
      void send(unsigned num, T* arr, unsigned recvRank, int tag);

/**
 * Receive an array from another rank.
 *
 * @param num number of elements to recieve.
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

    private:
/** MPI communicator */
      MPI_Comm _comm;
/** Flag to indicate if we are sending message */
      bool _sending;
/** Request object */
      MPI_Request _sendRequest;
/** Operators for allReduce */
      Lucee::MpiMsgOp _ops;
  };

}

#endif // LC_PAR_MPI_TMPL_H
