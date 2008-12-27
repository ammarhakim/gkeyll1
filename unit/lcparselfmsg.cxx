/**
 * @file	lcparselfmsg.cxx
 *
 * @brief	Unit tests for Lucee::ParSelfMsg class.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#include <lcparselfmsg.h>
#include <lctest.h>

#include <iostream>
#include <vector>

using namespace std;

void
test_a(Lucee::ParSelfMsg& selfMsg)
{
  LC_ASSERT("Testing number of processes", 
    selfMsg.numProcs() == 1);
  LC_ASSERT("Testing if rank is correct",
    selfMsg.rank() == 0);
}

void
test_all_reduce(Lucee::ParSelfMsg& selfMsg)
{
  unsigned val = 2 + selfMsg.rank();
  unsigned minVal, maxVal;
  selfMsg.allReduce(1, &val, &minVal, Lucee::PAR_MSG_MIN);
  LC_ASSERT("Testing if allReduce min worked",
    minVal == 2);
}

int
main(int argc, char *argv[])
{
  Lucee::ParSelfMsg selfMsg;
    
  LC_BEGIN_TESTS("lcparselfmsg");
  test_a(selfMsg);
  test_all_reduce(selfMsg);
  LC_END_TESTS;
}
