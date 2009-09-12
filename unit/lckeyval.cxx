/**
 * @file	lckeyval.cxx
 *
 * @brief	Unit tests for Lucee::KeyVal class
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// lucee includes
#include <LcKeyVal.h>
#include <LcTest.h>

void
test_1()
{
  Lucee::KeyVal kv;
  kv.add("nx", 10);
  kv.add("length", 1.0);

  std::vector<int> cells(2);
  cells[0] = 10; cells[1] = 20;
  kv.addVec("cells", cells);

  LC_ASSERT("Testing if fetching value works", kv.get<int>("nx") == 10);
  LC_ASSERT("Testing if fetching value works", kv.get<double>("length") == 1.0);

  LC_RAISES("Testing if exception is throw when getting non-exisitent data",
    kv.get<int>("length"), Lucee::Except);

  std::vector<int> gotCells = kv.getVec<int>("cells");
  LC_ASSERT("Testing length of vector", gotCells.size() == cells.size());
  for (unsigned i=0; i<cells.size(); ++i)
    LC_ASSERT("Testing if vector gotten properly", gotCells[i] == cells[i]);
}

int
main(void) 
{
  LC_BEGIN_TESTS("lckeyval");
  test_1();
  LC_END_TESTS;
}
