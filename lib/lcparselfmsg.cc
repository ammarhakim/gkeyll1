/**
 * @file	lcparselfmsg.cc
 *
 * @brief	Class for messaging in serial.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

// lib include
#include <lcparselfmsg.h>

namespace Lucee
{
  ParSelfMsg::ParSelfMsg()
  {
    // add templated messengers to the base class
    this->addMsg<char>
      ( new ParSelfTmpl<char>() );
    this->addMsg<unsigned char>
      ( new ParSelfTmpl<unsigned char>() );
    this->addMsg<short>
      ( new ParSelfTmpl<short>() );
    this->addMsg<unsigned short>
      ( new ParSelfTmpl<unsigned short>() );
    this->addMsg<int>
      ( new ParSelfTmpl<int>() );
    this->addMsg<unsigned int>
      ( new ParSelfTmpl<unsigned int>() );
    this->addMsg<long>
      ( new ParSelfTmpl<long>() );
    this->addMsg<unsigned long>
      ( new ParSelfTmpl<unsigned long>() );
    this->addMsg<float>
      ( new ParSelfTmpl<float>() );
    this->addMsg<double>
      ( new ParSelfTmpl<double>() );
    this->addMsg<long double>
      ( new ParSelfTmpl<long double>() );
    this->addMsg<long long int>
      ( new ParSelfTmpl<long long int>() );
  }

  ParMsgBase*
  ParSelfMsg::createSubComm(const std::vector<int>& ranks) 
  {
    return new ParSelfMsg();
  }

  void * 
  ParSelfMsg::finishRecv(MsgStatus lcms)
  {
    return 0;
  }

  bool
  ParSelfMsg::checkRecv(MsgStatus lcms)
  {
    return true;
  }
}
