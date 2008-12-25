/**
 * @file	lcparmpitmpl.cc
 *
 * @brief	Class for parallel messaging using MPI library.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#include <lcdatatypes.h>
#include <lcparmpitmpl.h>

namespace Lucee
{

  template<typename T>
  void 
  ParMpiTmpl<T>::send(const std::vector<T>& arr, unsigned recvRank, int tag) 
  {
    unsigned num = arr.size();
    this->resizeSendBuff(num);
    for (unsigned i=0; i<num; ++i)
      this->_sendBuff[i] = arr[i];

    MPI_Send( (void*) this->_sendBuff,
      num,
      MpiTraits<T>::mpiType(),
      recvRank,
      tag,
      _comm );
  }

  template<typename T>
  void 
  ParMpiTmpl<T>::send(unsigned num, T* arr, unsigned recvRank, int tag) 
  {
    MPI_Send( (void*) arr,
      num,
      MpiTraits<T>::mpiType(),
      recvRank,
      tag,
      _comm );
  }

  template<typename T>
  void 
  ParMpiTmpl<T>::recv(unsigned num, std::vector<T>& array, unsigned sendRank, int tag) 
  {
    MPI_Status status;
    this->resizeRecvBuff(num);
    // do a blocking receive
    MPI_Recv( (void*) this->_recvBuff,
      num,
      MpiTraits<T>::mpiType(),
      sendRank,
      tag,
      _comm,
      &status );
    // stick stuff into the vector
    for (unsigned i=0; i<num; ++i)
      array.push_back(this->_recvBuff[i]);
  }

  template<typename T>
  void 
  ParMpiTmpl<T>::recv(unsigned num, T* array, unsigned sendRank, int tag) 
  {
    MPI_Status status;
    // do a blocking receive
    MPI_Recv( (void*) array,
      num,
      MpiTraits<T>::mpiType(),
      sendRank,
      tag,
      _comm,
      &status );
  }

  template<typename T>
  MsgStatus 
  ParMpiTmpl<T>::startRecv(unsigned num, unsigned sendRank, int tag)
  {
    MpiMsgStatus_v *ms = new MpiMsgStatus_v();
    ms->data = new T[num];
    MPI_Irecv(ms->data,
      num,
      MpiTraits<T>::mpiType(),
      sendRank,
      tag,
      _comm,
      &ms->request);
    return ms;
  }

  template<typename T>
  void
  ParMpiTmpl<T>::allReduce(unsigned num, T* sendBuff, T* recvBuff, Lucee::ParMsgOpCode op)
  {
    MPI_Allreduce(sendBuff,
      recvBuff,
      num,
      MpiTraits<T>::mpiType(),
      _ops.getOp(op),
      _comm);
  }

// instantiate classes for all types define in wxdatatypes.h
  template class ParMpiTmpl<char>;
  template class ParMpiTmpl<unsigned char>;
  template class ParMpiTmpl<short>;
  template class ParMpiTmpl<unsigned short>;
  template class ParMpiTmpl<int>;
  template class ParMpiTmpl<unsigned int>;
  template class ParMpiTmpl<long>;
  template class ParMpiTmpl<unsigned long>;
  template class ParMpiTmpl<float>;
  template class ParMpiTmpl<double>;
  template class ParMpiTmpl<long double>;
  template class ParMpiTmpl<long long int>;

}
