/**
 * @file	lccurleval.cxx
 *
 * @brief	Unit tests for Lucee::CurlEval class
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// lucee includes
#include <LcCurlEval.h>
#include <LcTest.h>

void
resetCv(double cv[3])
{
  cv[0] = 2;
  cv[1] = 3;
  cv[2] = 4;
}

void
test_1()
{
  double delta = 1.0;
  double v[3] = {1, 2, 3};
  double vr[3] = {11, 12, 13};
  double cv[3];

  resetCv(cv);
// x-direction update
  Lucee::CurlEval<0>::eval(delta, v, vr, cv);
  LC_ASSERT("Testing curl calculation", cv[0] == 2.0);
  LC_ASSERT("Testing curl calculation", cv[1] == 3.0-10.0);
  LC_ASSERT("Testing curl calculation", cv[2] == 4.0+10.0);

  resetCv(cv);
// y-direction update
  Lucee::CurlEval<1>::eval(delta, v, vr, cv);
  LC_ASSERT("Testing curl calculation", cv[0] == 2.0+10);
  LC_ASSERT("Testing curl calculation", cv[1] == 3.0);
  LC_ASSERT("Testing curl calculation", cv[2] == 4.0-10.0);

  resetCv(cv);
// z-direction update
  Lucee::CurlEval<2>::eval(delta, v, vr, cv);
  LC_ASSERT("Testing curl calculation", cv[0] == 2.0-10);
  LC_ASSERT("Testing curl calculation", cv[1] == 3.0+10);
  LC_ASSERT("Testing curl calculation", cv[2] == 4.0);
}

int
main(void)
{
  LC_BEGIN_TESTS("lccurleval");
  test_1();
  LC_END_TESTS;
}
