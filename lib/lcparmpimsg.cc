/**
 * @file	lcparmpimsg.cc
 *
 * @brief	Class for parallel messaging using MPI library.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

// lib includes
#include <lcparmpimsg.h>

namespace Lucee
{

  ParMpiMsg::ParMpiMsg()
    : ParMsgBase(0, MPI_ANY_TAG), _comm(MPI_COMM_WORLD) 
  {
    // add templated messengers to the base class
    this->addMsg<char>
      ( new ParMpiTmpl<char>(_comm) );
    this->addMsg<unsigned char>
      ( new ParMpiTmpl<unsigned char>(_comm) );
    this->addMsg<short>
      ( new ParMpiTmpl<short>(_comm) );
    this->addMsg<unsigned short>
      ( new ParMpiTmpl<unsigned short>(_comm) );
    this->addMsg<int>
      ( new ParMpiTmpl<int>(_comm) );
    this->addMsg<unsigned int>
      ( new ParMpiTmpl<unsigned int>(_comm) );
    this->addMsg<long>
      ( new ParMpiTmpl<long>(_comm) );
    this->addMsg<unsigned long>
      ( new ParMpiTmpl<unsigned long>(_comm) );
    this->addMsg<float>
      ( new ParMpiTmpl<float>(_comm) );
    this->addMsg<double>
      ( new ParMpiTmpl<double>(_comm) );
    this->addMsg<long double>
      ( new ParMpiTmpl<long double>(_comm) );
    this->addMsg<long long int>
      ( new ParMpiTmpl<long long int>(_comm) );
  }

  ParMpiMsg::ParMpiMsg(ParMpiMsg *parent, MPI_Comm comm)
    : ParMsgBase(0, MPI_ANY_TAG, parent) {
    // add templated messengers to the base class
    this->addMsg<char>
      ( new ParMpiTmpl<char>(_comm) );
    this->addMsg<unsigned char>
      ( new ParMpiTmpl<unsigned char>(_comm) );
    this->addMsg<short>
      ( new ParMpiTmpl<short>(_comm) );
    this->addMsg<unsigned short>
      ( new ParMpiTmpl<unsigned short>(_comm) );
    this->addMsg<int>
      ( new ParMpiTmpl<int>(_comm) );
    this->addMsg<unsigned int>
      ( new ParMpiTmpl<unsigned int>(_comm) );
    this->addMsg<long>
      ( new ParMpiTmpl<long>(_comm) );
    this->addMsg<unsigned long>
      ( new ParMpiTmpl<unsigned long>(_comm) );
    this->addMsg<float>
      ( new ParMpiTmpl<float>(_comm) );
    this->addMsg<double>
      ( new ParMpiTmpl<double>(_comm) );
    this->addMsg<long double>
      ( new ParMpiTmpl<long double>(_comm) );
    this->addMsg<long long int>
      ( new ParMpiTmpl<long long int>(_comm) );
  }

  ParMsgBase*
  ParMpiMsg::createSubComm(const std::vector<int>& ranks) 
  {
    MPI_Group old_group, new_group;
    MPI_Comm new_comm;
    int *my_ranks = new int[ranks.size()];
    for (unsigned i=0; i<ranks.size(); ++i)
      my_ranks[i] = ranks[i];
    MPI_Comm_group(_comm, &old_group);
    MPI_Group_incl(old_group, ranks.size(), my_ranks, &new_group);
    MPI_Comm_create(_comm, new_group, &new_comm);
    delete [] my_ranks;

    if (new_comm == MPI_COMM_NULL)
      return 0;
    ParMsgBase *c = new ParMpiMsg(this, new_comm);
    return c;
  }

  void * 
  ParMpiMsg::finishRecv(MsgStatus wxms)
  {
    MpiMsgStatus_v *ms = static_cast<MpiMsgStatus_v *>(wxms);
    MPI_Status status;
    MPI_Wait(&ms->request, &status);
    void *data = ms->data;
    delete wxms;
    return data;
  }

  bool
  ParMpiMsg::checkRecv(MsgStatus wxms)
  {
    MpiMsgStatus_v *ms = static_cast<MpiMsgStatus_v *>(wxms);
    MPI_Status status;
    int flag;
    MPI_Test(&ms->request, &flag, &status);
    return (flag != 0) ? true : false;
  }

}
