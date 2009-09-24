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

void test_kv(const Lucee::KeyVal& kv)
{
  std::vector<int> cells(2);
  cells[0] = 10; cells[1] = 20;

  LC_ASSERT("Testing number of integers set", kv.getNum<int>() == 3);
  LC_ASSERT("Testing number of doubles set", kv.getNum<double>() == 3);
  LC_ASSERT("Testing number of vector of ints set", 
    kv.getNum<std::vector<int> >() == 1);

  LC_ASSERT("Testing if fetching value works", kv.get<int>("nx") == 10);
  LC_ASSERT("Testing if fetching value works", kv.get<int>("ny") == 20);
  LC_ASSERT("Testing if fetching value works", kv.get<int>("nz") == 30);
  LC_ASSERT("Testing if fetching value works", kv.get<double>("lengthx") == 1.0);
  LC_ASSERT("Testing if fetching value works", kv.get<double>("lengthy") == 2.0);
  LC_ASSERT("Testing if fetching value works", kv.get<double>("lengthz") == 3.0);

  std::vector<int> gotCells = kv.get<std::vector<int> >("cells");
  LC_ASSERT("Testing length of vector", gotCells.size() == cells.size());
  for (unsigned i=0; i<cells.size(); ++i)
    LC_ASSERT("Testing if vector gotten properly", gotCells[i] == cells[i]);

  LC_ASSERT("Testing to see if has() works", kv.has<int>("nx") == true);
  LC_ASSERT("Testing to see if has() works", kv.has<int>("ny") == true);
  LC_ASSERT("Testing to see if has() works", kv.has<int>("nz") == true);
  LC_ASSERT("Testing to see if has() works", kv.has<double>("lengthx") == true);
  LC_ASSERT("Testing to see if has() works", kv.has<double>("lengthy") == true);
  LC_ASSERT("Testing to see if has() works", kv.has<double>("lengthz") == true);
  LC_ASSERT("Testing to see if has() works", kv.has<std::vector<int> >("cells") == true);
  LC_ASSERT("Testing to see if has() works", kv.has<double>("foo") == false);

  LC_RAISES("Testing if exception is throw when getting non-existent data",
    kv.get<int>("length"), Lucee::Except);

// check iterators
  kv.setToFirst<int>();
  for (unsigned i=0; i<kv.getNum<int>(); ++i)
  {
    std::pair<std::string, int> kvp = kv.getAndBump<int>();
// test if correct pair was returned
    LC_ASSERT("Testing if returned pair was correct", kv.has<int>(kvp.first));
    LC_ASSERT("Testing if returned pair was correct", kv.get<int>(kvp.first) == kvp.second);
  }

  kv.setToFirst<double>();
  for (unsigned i=0; i<kv.getNum<double>(); ++i)
  {
    std::pair<std::string, double> kvp = kv.getAndBump<double>();
// test if correct pair was returned
    LC_ASSERT("Testing if returned pair was correct", kv.has<double>(kvp.first));
    LC_ASSERT("Testing if returned pair was correct", kv.get<double>(kvp.first) == kvp.second);
  }

  kv.setToFirst<std::vector<int> >();
  for (unsigned i=0; i<kv.getNum<std::vector<int> >(); ++i)
  {
    std::pair<std::string, std::vector<int> > kvp = kv.getAndBump<std::vector<int> >();
// test if correct pair was returned
    LC_ASSERT("Testing if returned pair was correct", 
      kv.has<std::vector<int> >(kvp.first));
    std::vector<int> vec = kv.get<std::vector<int> >(kvp.first);
    for (unsigned k=0; k<vec.size(); ++k)
      LC_ASSERT("Testing if returned pair was correct", vec[k] == kvp.second[k]);
  }
}

void
test_1()
{
  Lucee::KeyVal kv;

// add stuff to set
  kv.add("nx", 10);
  kv.add("ny", 20);
  kv.add("nz", 30);
  kv.add("lengthx", 1.0);
  kv.add("lengthy", 2.0);
  kv.add("lengthz", 3.0);
  std::vector<int> cells(2);
  cells[0] = 10; cells[1] = 20;
  kv.add("cells", cells);

// check stuff
  test_kv(kv);

// make a copy
  Lucee::KeyVal kvc(kv);
  test_kv(kvc);

// make an assignment
  Lucee::KeyVal kva;
  kva = kv;
  test_kv(kva);

// make a duplicate
  Lucee::KeyVal kvd = kv.duplicate();
  test_kv(kvd);

// modify to see shallow copies were made
  kva.add("components", 5);
  LC_ASSERT("Testing if copy is shallow", kv.get<int>("components") == 5);
  LC_ASSERT("Testing if copy is shallow", kvc.get<int>("components") == 5);

// make sure duplicate copy did not change
  LC_RAISES("Ensuring that duplicate was not changed", kvd.get<int>("components"),
    Lucee::Except);
}

int
main(void) 
{
  LC_BEGIN_TESTS("lckeyval");
  test_1();
  LC_END_TESTS;
}
