/**
 * @file	lcparmpimsg.cxx
 *
 * @brief	Unit tests for Lucee::ParMpiMsg class.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#include <lcparmpimsg.h>
#include <lctest.h>

#include <iostream>
#include <vector>

using namespace std;

void
test_a(Lucee::ParMpiMsg& mpiMsg)
{
  if (mpiMsg.rank() == 0)
  {
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    LC_ASSERT("Testing number of processes", 
      mpiMsg.numProcs() == (unsigned) nproc);
  }
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  LC_ASSERT("Testing if rank is correct",
    mpiMsg.rank() == rank);

  // wait for everyone to catch up
  mpiMsg.barrier();
}

template<typename T>
void
test_untagged_send_recv(Lucee::ParMpiMsg& mpiMsg)
{
  if (mpiMsg.rank() != 0)
  {
    T sdata_1[10];
    for (unsigned i=0; i<10; ++i)
      sdata_1[i] = (T) mpiMsg.rank() + (T) i + 1;

    // send stuff to process 0 without tags
    mpiMsg.send(10, sdata_1, 0);
  }
  else
  {
    T rdata[100];
    for (unsigned i=1; i<mpiMsg.numProcs(); ++i)
    {
      // construct arrays to compare with
      T sdata_1[10];
      for (unsigned j=0; j<10; ++j)
        sdata_1[j] = (T) i + (T) j + 1;

      mpiMsg.recv(10, rdata, i);
      LC_ASSERT("Testing to see if message without tag was recieved properly",
        arraycmp(rdata, sdata_1, 10));
    }
  }

  // wait for everyone to catch up
  mpiMsg.barrier();
}

template<typename T>
void
test_tagged_send_recv(Lucee::ParMpiMsg& mpiMsg)
{
  if (mpiMsg.rank() != 0)
  {
    T sdata_1[10];
    for (unsigned i=0; i<10; ++i)
      sdata_1[i] = (T) mpiMsg.rank() + (T) i + 1;

    T sdata_2[10];
    for (unsigned i=0; i<10; ++i)
      sdata_2[i] = (T) mpiMsg.rank() + (T) i + 2;

    // send stuff to process 0
    mpiMsg.send(10, sdata_1, 0, 9);
    mpiMsg.send(10, sdata_2, 0, 10);
  }
  else
  {
    T rdata[100];
    for (unsigned i=1; i<mpiMsg.numProcs(); ++i)
    {
      // construct arrays to compare with
      T sdata_1[10];
      for (unsigned j=0; j<10; ++j)
        sdata_1[j] = (T) i + (T) j + 1;
      T sdata_2[10];
      for (unsigned j=0; j<10; ++j)
        sdata_2[j] = (T) i + (T) j + 2;

      mpiMsg.recv(10, rdata, i, 10);
      LC_ASSERT("Testing to see if message with tag 10 was recieved properly",
        arraycmp(rdata, sdata_2, 10));
            
      mpiMsg.recv(10, rdata, i, 9);
      LC_ASSERT("Testing to see if message with tag 9 was recieved properly",
        arraycmp(rdata, sdata_1, 10));
    }
  }

  // wait for everyone to catch up
  mpiMsg.barrier();
}

template<typename T>
void
test_non_block(Lucee::ParMpiMsg& mpiMsg)
{

  if (mpiMsg.rank() != 0)
  {
    T sdata_1[10];
    for (unsigned i=0; i<10; ++i)
      sdata_1[i] = (T) mpiMsg.rank() + (T) i + 1;

    T sdata_2[10];
    for (unsigned i=0; i<10; ++i)
      sdata_2[i] = (T) mpiMsg.rank() + (T) i + 2;

    // send stuff to process 0
    mpiMsg.send(10, sdata_1, 0, 9);
    mpiMsg.send(10, sdata_2, 0, 10);
  }
  else
  {
    for (unsigned i=1; i<mpiMsg.numProcs(); ++i)
    {
      // construct arrays to compare with
      T sdata_1[10];
      for (unsigned j=0; j<10; ++j)
        sdata_1[j] = (T) i + (T) j + 1;
      T sdata_2[10];
      for (unsigned j=0; j<10; ++j)
        sdata_2[j] = (T) i + (T) j + 2;

      T *rdata;

      Lucee::MsgStatus ms10 = mpiMsg.template
        startRecv<T>(10, i, 10);
      rdata = (T*) mpiMsg.finishRecv(ms10);
      LC_ASSERT("Testing to see if message with tag 10 was non-blocking recieved properly",
        arraycmp(rdata, sdata_2, 10));
      LC_ASSERT("Making sure that checkRecv returns true",
        mpiMsg.checkRecv(ms10) == true);
      delete [] rdata;

      Lucee::MsgStatus ms9 = mpiMsg.template
        startRecv<T>(10, i, 9);
      rdata = (T*) mpiMsg.finishRecv(ms9);
      LC_ASSERT("Testing to see if message with tag 9 was non-blocking recieved properly",
        arraycmp(rdata, sdata_1, 10));
      LC_ASSERT("Making sure that checkRecv returns true",
        mpiMsg.checkRecv(ms9) == true);
      delete [] rdata;
    }
  }
  // wait for everyone to catch up
  mpiMsg.barrier();
}

void
test_all_reduce(Lucee::ParMpiMsg& mpiMsg)
{
  unsigned val = 2 + mpiMsg.rank();
  unsigned minVal, maxVal;
  mpiMsg.allReduce(1, &val, &minVal, Lucee::PAR_MSG_MIN);
  LC_ASSERT("Testing if allReduce min worked",
    minVal == 2);
}

int
main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);

  Lucee::ParMpiMsg mpiMsg;
    
  LC_MPI_BEGIN_TESTS("lcparmpimsg");

  test_a(mpiMsg);

  test_tagged_send_recv<char>(mpiMsg);
  test_tagged_send_recv<unsigned char>(mpiMsg);
  test_tagged_send_recv<short>(mpiMsg);
  test_tagged_send_recv<unsigned short>(mpiMsg);
  test_tagged_send_recv<int>(mpiMsg);
  test_tagged_send_recv<unsigned int>(mpiMsg);
  test_tagged_send_recv<long>(mpiMsg);
  test_tagged_send_recv<unsigned long>(mpiMsg);
  test_tagged_send_recv<float>(mpiMsg);
  test_tagged_send_recv<double>(mpiMsg);
  test_tagged_send_recv<long double>(mpiMsg);
  test_tagged_send_recv<long long int>(mpiMsg);

  test_untagged_send_recv<char>(mpiMsg);
  test_untagged_send_recv<unsigned char>(mpiMsg);
  test_untagged_send_recv<short>(mpiMsg);
  test_untagged_send_recv<unsigned short>(mpiMsg);
  test_untagged_send_recv<int>(mpiMsg);
  test_untagged_send_recv<unsigned int>(mpiMsg);
  test_untagged_send_recv<long>(mpiMsg);
  test_untagged_send_recv<unsigned long>(mpiMsg);
  test_untagged_send_recv<float>(mpiMsg);
  test_untagged_send_recv<double>(mpiMsg);
  test_untagged_send_recv<long double>(mpiMsg);
  test_untagged_send_recv<long long int>(mpiMsg);

  test_non_block<char>(mpiMsg);
  test_non_block<unsigned char>(mpiMsg);
  test_non_block<short>(mpiMsg);
  test_non_block<unsigned short>(mpiMsg);
  test_non_block<int>(mpiMsg);
  test_non_block<unsigned int>(mpiMsg);
  test_non_block<long>(mpiMsg);
  test_non_block<unsigned long>(mpiMsg);
  test_non_block<float>(mpiMsg);
  test_non_block<double>(mpiMsg);
  test_non_block<long double>(mpiMsg);
  test_non_block<long long int>(mpiMsg);

  test_all_reduce(mpiMsg);

  LC_MPI_END_TESTS;

  MPI_Finalize();
}
