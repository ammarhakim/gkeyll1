/**
 * @file	lckeyvaltreelexer.cxx
 *
 * @brief	Unit tests for Lucee::KeyValueTreeLexer class
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#include <lckeyvaltreelexer.h>
#include <lctest.h>

// std include
#include <sstream>

using namespace std;

void
test_a()
{
  std::istringstream iss("<lucee> </lucee>");
  Lucee::KeyValTreeLexer<float> lexer(iss);

  LC_ASSERT("Testing token in input stream",
    lexer.YYLex() == Lucee::LC_LEFT_ANGLE);
 
  LC_ASSERT("Testing token in input stream",
    lexer.YYLex() == Lucee::LC_ID);
  LC_ASSERT("Testing token string in input stream",
    lexer.YYText() == "lucee");

  LC_ASSERT("Testing token in input stream",
    lexer.YYLex() == Lucee::LC_RIGHT_ANGLE);

  LC_ASSERT("Testing token in input stream",
    lexer.YYLex() == Lucee::LC_LEFT_ANGLE_FRONT_SLASH);

  LC_ASSERT("Testing token in input stream",
    lexer.YYLex() == Lucee::LC_ID);
  LC_ASSERT("Testing token string in input stream",
    lexer.YYText() == "lucee");

  LC_ASSERT("Testing token in input stream",
    lexer.YYLex() == Lucee::LC_RIGHT_ANGLE);
}

int
main(void)
{
  LC_BEGIN_TESTS("lckeyvaltreelexer");
  test_a();
  LC_END_TESTS;
}
