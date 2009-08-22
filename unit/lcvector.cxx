/**
 * @file	lcvector.cxx
 *
 * @brief	Unit tests for Lucee::Vector class
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// lucee includes
#include <LcVector.h>
#include <LcTest.h>

// std includes
#include <cmath>

void
test_1()
{
  Lucee::Vector<double> vec(10);

  LC_ASSERT("Testing size of vector", vec.getLength() == 10);
  vec = 10.5;

  for (unsigned i=0; i<vec.getLength(); ++i)
    LC_ASSERT("Testing value of fixed vector", vec[i] == 10.5);

  for (unsigned i=0; i<vec.getLength(); ++i)
    vec[i] = i;

  for (unsigned i=0; i<vec.getLength(); ++i)
    LC_ASSERT("Testing value of fixed vector", vec[i] == i);
}

int
main(void)
{
  LC_BEGIN_TESTS("lcvector");
  test_1();
  LC_END_TESTS;
}
