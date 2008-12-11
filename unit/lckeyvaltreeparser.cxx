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

int
main(void)
{
  LC_BEGIN_TESTS("lckeyvaltreeparser");
  test_a();
  LC_END_TESTS;
}
