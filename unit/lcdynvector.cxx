/**
 * @file	lcdynvector.cxx
 *
 * @brief	Unit tests for Lucee::Array class
 */

// lucee includes
#include <LcDynVector.h>
#include <LcTest.h>

void
test_1()
{
  Lucee::DynVector<double> dynVec(3);
  std::vector<double> xp(3);
  xp[0] = 1.0; xp[1] = 2.0; xp[2] = 3.0;
  dynVec.appendData(0.1, xp);

  std::vector<double> val = dynVec.getLastInsertedData();
  LC_ASSERT("Testing if last inserted data was fetched", val[0] == 1.0);
  LC_ASSERT("Testing if last inserted data was fetched", val[1] == 2.0);
  LC_ASSERT("Testing if last inserted data was fetched", val[2] == 3.0);
  LC_ASSERT("Testing if time was fetched", dynVec.getLastInsertedTime() == 0.1);

  xp[0] = 1.0; xp[1] = 2.0; xp[2] = 3.0;
  dynVec.appendData(0.2, xp);

  val = dynVec.getLastInsertedData();
  LC_ASSERT("Testing if last inserted data was fetched", val[0] == 1.0);
  LC_ASSERT("Testing if last inserted data was fetched", val[1] == 2.0);
  LC_ASSERT("Testing if last inserted data was fetched", val[2] == 3.0);
  LC_ASSERT("Testing if time was fetched", dynVec.getLastInsertedTime() == 0.2);

  LC_ASSERT("Testing if numComponents is correct",
    dynVec.getNumComponents() == 3);
  LC_ASSERT("Testing if size is correct", dynVec.getSize() == 2);

  dynVec.removeLastInsertedData();
  val = dynVec.getLastInsertedData();
  LC_ASSERT("Testing if last inserted data was fetched", val[0] == 1.0);
  LC_ASSERT("Testing if last inserted data was fetched", val[1] == 2.0);
  LC_ASSERT("Testing if last inserted data was fetched", val[2] == 3.0);
  LC_ASSERT("Testing if time was fetched", dynVec.getLastInsertedTime() == 0.1);

  std::vector<double> wrong(2);
  LC_RAISES("Testing to see if incorrect data can be appended",
    dynVec.appendData(0.1, wrong), Lucee::Except);

  dynVec.clear();
  LC_ASSERT("Testing if size is correct", dynVec.getSize() == 0);
}

int
main(int argc, char **argv) 
{
  LC_BEGIN_TESTS("lcdynvector");
  test_1();
  LC_END_TESTS;
}
