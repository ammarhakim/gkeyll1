/**
 * @file	lcexprparser.cxx
 *
 * @brief	Unit tests for Lucee::ExprParser class
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

// lib includes
#include <lctest.h>
#include <lcexprparser.h>

// std includes
#include <vector>

void
test_a()
{
  Lucee::ExprParser ep;
  ep.appendIndVar("x");
  ep.appendIndVar("y");
  ep.addConstant("c0", 1.0);
  unsigned res1 = ep.addExpr("x+y+c0");
  ep.setup();

  std::vector<double> iv;
  iv.push_back(1.0); // x
  iv.push_back(2.0); // y
  
  ep.eval(iv);
  LC_ASSERT("Testing if expression result is correct",
    ep.result(res1) == 4.0);
}

int
main(void)
{
  LC_BEGIN_TESTS("lcexprparser");
  test_a();
  LC_END_TESTS;
}
