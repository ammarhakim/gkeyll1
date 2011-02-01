/**
 * @file	lcluahelp.cxx
 *
 * @brief	Unit tests for Luahelp and ConstLuahelp classes
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// lucee includes
#include <LcLuaHelp.h>
#include <LcTest.h>

// Loki includes
#include <loki/Singleton.h>

// std includes
#include <string>

class Base
{
  public:
    static std::string describe()
    {
      return "Base class help\n";
    }
};

class Derived : public Base
{
  public:
    static std::string describe()
    {
      return "Derived class help\n";
    }
};

void
test_1()
{
}

int
main(void) 
{
  LC_BEGIN_TESTS("lcluahelp");
  test_1();
  LC_END_TESTS;
}

