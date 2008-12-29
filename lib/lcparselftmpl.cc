/**
 * @file	lcparselftmpl.cc
 *
 * @brief	Class for messaging in serial.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

// lib includes
#include <lcdatatypes.h>
#include <lcparselftmpl.h>

namespace Lucee
{
  template<typename T>
  void 
  ParSelfTmpl<T>::send(const std::vector<T>& arr, unsigned recvRank, int tag) 
  {
  }

  template<typename T>
  void 
  ParSelfTmpl<T>::send(unsigned num, T* arr, unsigned recvRank, int tag) 
  {
  }

  template<typename T>
  void 
  ParSelfTmpl<T>::recv(unsigned num, std::vector<T>& array, unsigned sendRank, int tag) 
  {
  }

  template<typename T>
  void 
  ParSelfTmpl<T>::recv(unsigned num, T* array, unsigned sendRank, int tag) 
  {
  }

  template<typename T>
  MsgStatus 
  ParSelfTmpl<T>::startRecv(unsigned num, unsigned sendRank, int tag)
  {
    return new SelfMsgStatus_v();
  }

  template<typename T>
  void
  ParSelfTmpl<T>::allReduce(unsigned num, T* sendBuff, T* recvBuff, Lucee::ParMsgOpCode op)
  {
    for (unsigned i=0; i<num; ++i)
      recvBuff[i] = sendBuff[i];
  }

// instantiate classes for all types define in lcdatatypes.h
  template class ParSelfTmpl<char>;
  template class ParSelfTmpl<unsigned char>;
  template class ParSelfTmpl<short>;
  template class ParSelfTmpl<unsigned short>;
  template class ParSelfTmpl<int>;
  template class ParSelfTmpl<unsigned int>;
  template class ParSelfTmpl<long>;
  template class ParSelfTmpl<unsigned long>;
  template class ParSelfTmpl<float>;
  template class ParSelfTmpl<double>;
  template class ParSelfTmpl<long double>;
  template class ParSelfTmpl<long long int>;
}
