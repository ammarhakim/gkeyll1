/**
 * @file	lckeyvaltreeparser.cxx
 *
 * @brief	Unit tests for Lucee::KeyValueTreeParser class
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#include <lckeyvaltreeparser.h>
#include <lckeyvaltree.h>
#include <lctest.h>

// std include
#include <sstream>

using namespace std;

void
test_a()
{
  std::istringstream iss("<lucee> val = 36 </lucee>");

  Lucee::KeyValTree kvt;
  Lucee::KeyValTreeParser<double> parser(iss, kvt);
  parser.parse();

  LC_ASSERT("Testing name of KeyValueSet",
    kvt.getName() == "lucee");
  LC_ASSERT("Testing stored value",
    kvt.get<int>("val") == 36);
}

void
test_b()
{
  std::istringstream iss("<lucee> vals = [1, 2] </lucee>");

  Lucee::KeyValTree kvt;
  Lucee::KeyValTreeParser<double> parser(iss, kvt);
  parser.parse();

  LC_ASSERT("Testing name of KeyValueSet",
    kvt.getName() == "lucee");
  std::vector<int> vals = kvt.get<std::vector<int> >("vals");
  LC_ASSERT("Testing if list size if proper",
    vals.size() == 2);
  LC_ASSERT("Testing values gotten from list",
    vals[0] == 1);
  LC_ASSERT("Testing values gotten from list",
    vals[1] == 2);
}

int
main(void)
{
  LC_BEGIN_TESTS("lckeyvaltreeparser");
  test_a();
  test_b();
  LC_END_TESTS;
}
