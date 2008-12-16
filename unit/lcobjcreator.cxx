/**
 * @file	lccreator.cxx
 *
 * @brief	Unit tests for Lucee::ObjCreator and related classes
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#include <lcobjregister.h>
#include <lcobjcreator.h>
#include <lctest.h>

class Base
{
  public:
    virtual ~Base() {}

    virtual std::string getMessage() = 0;
};

class Derived : public Base 
{
  public:
    std::string getMessage()
    {
      return "I am Dervied";
    }
};

Lucee::ObjRegister<Derived, Base> base("base");

void
test_a()
{
  Base *d = Lucee::ObjCreator<Base>::getNew("base");
  LC_ASSERT("Testing if we got proper type",
    d->getMessage() == "I am Dervied");
}

int
main(void)
{
  LC_BEGIN_TESTS("lcobjcreator");
  test_a();
  LC_END_TESTS;
}

