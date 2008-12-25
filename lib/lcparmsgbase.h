/**
 * @file	lcparmsgbase.h
 *
 * @brief	Base class for parallel message passing classes.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_PAR_MSG_BASE_H
#define LC_PAR_MSG_BASE_H

// lib includes
#include <lcparmsgtmplbase.h>
#include <lcdatatypes.h>
#include <lcexcept.h>

// etc includes
#include <loki/HierarchyGenerators.h>

// std includes
#include <vector>

namespace Lucee
{

/**
 * Provides an abstract interface for message based communication
 * between different processes.
 */
  class ParMsgBase 
  {
    public:
/**
 * Destructor
 */
      virtual ~ParMsgBase() 
      {
      }

/**
 * Rank of process
 * 
 * @return this rank
 */
      virtual int rank() const = 0;

/**
 * Number of processes taking part in messaging
 *
 * @return num of processes
 */
      virtual unsigned numProcs() const = 0;

/**
 * Get parent communicating processor group
 *
 * @return parent proccesor group
 */
      ParMsgBase* parent() const 
      {
        return _parent;
      }

/**
 * Split communicator into a child communicator
 *
 * @param ranks list of processors in old communicator which are to be
 * in the new communicator
 * @return new communicator
 */
      virtual ParMsgBase* createSubComm(const std::vector<int>& ranks) = 0;


/**
 * Block till all processes hit this barrier
 */
      virtual void barrier() const = 0;

/**
 * Send a std::vector to another rank.
 * 
 * @param array std::vector of data being sent.
 * @param recvRank rank that will receive array.
 * @param tag Tag to attache to this message.
 */
      template<typename T>
      void send(const std::vector<T>& array, unsigned recvRank, int tag=-1) 
      {
        this->_getMsg<T>()->send(array, recvRank, (tag == -1) ? _sendTag : tag);
      }

/**
 * Send an array to another rank.
 * 
 * @param num number of elements to send.
 * @param array array of length 'num' of data being sent.
 * @param recvRank rank that will receive array.
 * @param tag Tag to attache to this message.
 */
      template<typename T>
      void send(unsigned num, T* array, unsigned recvRank, int tag=-1) 
      {
        this->_getMsg<T>()->send(num, array, recvRank, (tag == -1) ? _sendTag : tag);
      }

/**
 * Receive a std::vector from another rank.
 * 
 * @param num number of elements to receive.
 * @param array array that is filled with received values.
 * @param sendRank rank that sent array.
 * @param tag Tag to attache to this message.
 */
      template<typename T>
      void recv(int num, std::vector<T>& array, unsigned sendRank, int tag=-1) 
      {
        this->_getMsg<T>()->recv(num, array, sendRank, (tag == -1) ? _recvTag : tag);
      }

/**
 * Receive an array from another rank.
 * 
 * @param num number of elements to receive.
 * @param array array that is filled with received values.
 * @param sendRank rank that sent array.
 * @param tag Tag to attache to this message.
 */
      template<typename T>
      void recv(unsigned num, T* array, unsigned sendRank, int tag=-1) 
      {
        this->_getMsg<T>()->recv(num, array, sendRank, (tag == -1) ? _recvTag : tag);
      }

/**
 * Receive an array from another rank. This is a non-blocking call.
 * 
 * @param num number of elements to receive.
 * @param sendRank rank we want to receive from.
 * @param tag Tag to attache to this message.
 * @return message status.
 */
      template<typename T>
      MsgStatus startRecv(unsigned num,  unsigned sendRank, int tag=-1) 
      {
        return this->_getMsg<T>()->startRecv(
            num, sendRank, (tag == -1) ? _recvTag : tag);
      }

/**
 * Finish the receive started by startRecv and return a pointer to the
 * data recieved. The calling function owns the pointer and is
 * resposible for freeing it.
 *
 * @param ms message status object returned by startRecv
 * @return pointer to the data recieved. Caller owns the returned pointer
 */
      virtual void * finishRecv(MsgStatus ms) = 0;

/**
 * Check status of recieve started by a startRecv. If this call
 * returns true then finishRecv can be called to immediately recieve
 * the data.
 *
 * @param ms message status object returned by startRecv
 * @return true it recieve has been completed, false otherwise
 */
      virtual bool checkRecv(MsgStatus ms) = 0;

/**
 * Reduce data to all ranks
 * 
 * @param num Number of elements being reduced
 * @param sendBuff buffer to reduce
 * @param recvBuff buffer to recieve reduced data
 * @param op operation to perform. These are one of those listed in ParMsgOps enum
 */
      template<typename T>
      void allReduce(unsigned num, T* sendBuff, T* recvBuff, Lucee::ParMsgOps op) 
      {
        this->_getMsg<T>()->allReduce(num, sendBuff, recvBuff, op);
      }

    protected:
/**
 * Constructor. This is protected so only children can make instances.
 */
      ParMsgBase(int sendTag=0, int recvTag=0, ParMsgBase* parent=0)
        : _sendTag(sendTag), _recvTag(recvTag), _parent(parent) 
      {
      }

/**
 * Add a new messager : the derived class should call this to setup
 * ParMsgBase properly
 */
      template <typename T>
      void addMsg(Lucee::ParMsgTmplBase<T> *b) 
      {
        Loki::Field<T>(_msgTypeMap)._msg = b;
      }

    private:

/** 
 * Copy constructor: private to prevent use.
 *
 * @param rhs Object to copy.
 */
      ParMsgBase(const ParMsgBase& rhs);

/** 
 * Assignment operator: private to prevent use .
 *
 * @param msg Object to assign from.
 * @return reference to this object.
 */
      ParMsgBase& operator=(const ParMsgBase& msg);

/** Tag for sent message */
      int _sendTag;
/** Tag for recv message */
      int _recvTag;
/** Pointer to parent messenger object */
      ParMsgBase *_parent;

/**
 * Get a messager object with the proper type
 */
      template <typename T>
      ParMsgTmplBase<T>* _getMsg() 
      {
        ParMsgTmplBase<T> *r = Loki::Field<T>(_msgTypeMap)._msg;
        if (r) return r;
        Lucee::Except ex;
        ex << "Message type not set properly"; 
        throw ex;
      }

    public:

/** Container class for all messagers */
      template<typename T>
      struct MsgContainer 
      {
/** Create a new message container */
          MsgContainer() 
            : _msg(0) 
          {
          }

/** Delete the message contsainer */
          virtual ~MsgContainer() 
          {
            delete _msg;
          }
/**  Pointer to messaging object: points to derived class of ParMsgTmplBase<T> */
          ParMsgTmplBase<T> *_msg;
      };

/** Type definition of map of types to message containers */
      typedef Loki::GenScatterHierarchy<DataTypes_t, MsgContainer> MsgTypeMap_t;

/** Map of types to message containers */
      MsgTypeMap_t _msgTypeMap;
  };
}

#endif // LC_PAR_MSG_BASE_H
