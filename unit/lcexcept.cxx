/**
 * @file	lcexcept.cxx
 *
 * @brief	Unit tests for Lucee::Except class
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#include <lcexcept.h>
#include <lctest.h>

#include <string.h> // for c-style strcmp

using namespace std;

void
test_a()
{
  Lucee::Except ex("This is an exception message");
  LC_ASSERT("Testing if exception message set properly",
    ex.str() == "This is an exception message");

  LC_ASSERT("Testing if exception message set properly",
    strcmp(ex.what(), "This is an exception message") == 0);

  Lucee::Except ex1(ex);
  LC_ASSERT("Testing if copy ctor works",
    ex1.str() == "This is an exception message");

  Lucee::Except ex2("Hello");
  ex2 = ex1;
  LC_ASSERT("Testing if copy ctor works",
    ex2.str() == "This is an exception message");

  ex << " with an addional clause";
  LC_ASSERT("Testing if exception message set properly",
    ex.str() == "This is an exception message with an addional clause");
  LC_ASSERT("Testing if exception message set properly",
    strcmp(ex.what(), "This is an exception message with an addional clause") == 0);
}

int
main(void)
{
  LC_BEGIN_TESTS("lcexcept");
  test_a();
  LC_END_TESTS;
}
