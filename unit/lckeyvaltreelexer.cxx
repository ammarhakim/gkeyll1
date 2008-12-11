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

void
test_b()
{
  std::istringstream is1("val = 22");
  Lucee::KeyValTreeLexer<float> lex1(is1);

  LC_ASSERT("Testing token in input stream",
    lex1.YYLex() == Lucee::LC_ID);
  LC_ASSERT("Teting token string in input stream",
    lex1.YYText() == "val");

  LC_ASSERT("Testing token in input stream",
    lex1.YYLex() == Lucee::LC_EQUAL);

  LC_ASSERT("Testing token in input stream",
    lex1.YYLex() == Lucee::LC_INT);
  LC_ASSERT("Testing token value in input stream",
    lex1.integer() == 22);

  std::istringstream is2("val = \"Hello\"");
  Lucee::KeyValTreeLexer<float> lex2(is2);

  LC_ASSERT("Testing token in input stream",
    lex2.YYLex() == Lucee::LC_ID);
  LC_ASSERT("Teting token string in input stream",
    lex2.YYText() == "val");

  LC_ASSERT("Testing token in input stream",
    lex2.YYLex() == Lucee::LC_EQUAL);

  LC_ASSERT("Testing token in input stream",
    lex2.YYLex() == Lucee::LC_STRING);
  LC_ASSERT("Testing token value in input stream",
    lex2.YYText() == "Hello");

  std::istringstream is3("val = 3.14");
  Lucee::KeyValTreeLexer<double> lex3(is3);

  LC_ASSERT("Testing token in input stream",
    lex3.YYLex() == Lucee::LC_ID);
  LC_ASSERT("Teting token string in input stream",
    lex3.YYText() == "val");

  LC_ASSERT("Testing token in input stream",
    lex3.YYLex() == Lucee::LC_EQUAL);

  LC_ASSERT("Testing token in input stream",
    lex3.YYLex() == Lucee::LC_REAL);
  LC_ASSERT("Testing token value in input stream",
    lex3.real() == 3.14);

  std::istringstream is4("val = [1, 2, 3]");
  Lucee::KeyValTreeLexer<double> lex4(is4);

  LC_ASSERT("Testing token in input stream",
    lex4.YYLex() == Lucee::LC_ID);
  LC_ASSERT("Teting token string in input stream",
    lex4.YYText() == "val");

  LC_ASSERT("Testing token in input stream",
    lex4.YYLex() == Lucee::LC_EQUAL);

  LC_ASSERT("Testing token in input stream",
    lex4.YYLex() == Lucee::LC_LEFT_BOX);

  LC_ASSERT("Testing token in input stream",
    lex4.YYLex() == Lucee::LC_INT);
  LC_ASSERT("Testing token value in input stream",
    lex4.integer() == 1);

  LC_ASSERT("Testing token in input stream",
    lex4.YYLex() == Lucee::LC_COMMA);

  LC_ASSERT("Testing token in input stream",
    lex4.YYLex() == Lucee::LC_INT);
  LC_ASSERT("Testing token value in input stream",
    lex4.integer() == 2);

  LC_ASSERT("Testing token in input stream",
    lex4.YYLex() == Lucee::LC_COMMA);

  LC_ASSERT("Testing token in input stream",
    lex4.YYLex() == Lucee::LC_INT);
  LC_ASSERT("Testing token value in input stream",
    lex4.integer() == 3);

  LC_ASSERT("Testing token in input stream",
    lex4.YYLex() == Lucee::LC_RIGHT_BOX);
}

int
main(void)
{
  LC_BEGIN_TESTS("lckeyvaltreelexer");
  test_a();
  test_b();
  LC_END_TESTS;
}
