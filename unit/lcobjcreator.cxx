/**
 * @file	lccreator.cxx
 *
 * @brief	Unit tests for Lucee::ObjCreator and related classes
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// lib includes
#include <LcObjRegistry.h>
#include <LcObjCreator.h>
#include <LcTest.h>

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
      return "I am Derived";
    }
};

class Derived2 : public Base 
{
  public:
    std::string getMessage()
    {
      return "I am Derived2";
    }
};

Lucee::ObjRegistry<Base, Derived> derived("derived");
Lucee::ObjRegistry<Base, Derived2> derived2("derived2");

void
test_1()
{
  Base *d = Lucee::ObjCreator<Base>::getNew("derived");
  LC_ASSERT("Testing if we got proper type",
    d->getMessage() == "I am Derived");

  Base *d2 = Lucee::ObjCreator<Base>::getNew("derived2");
  LC_ASSERT("Testing if we got proper type",
    d2->getMessage() == "I am Derived2");

  std::vector<std::string> names
    = Lucee::ObjCreator<Base>::registeredNames();
  LC_ASSERT("Ensure we have two creators registered",
    names.size() == 2);

  LC_ASSERT("Testing if name is correct",
    names[0] == "derived");
  LC_ASSERT("Testing if name is correct",
    names[1] == "derived2");

  // test unregister method
  Lucee::ObjCreator<Base>::unregister("derived");
  LC_RAISES("Testing if exception is raised",
    Lucee::ObjCreator<Base>::getNew("derived"), std::exception);
}

int
main(void)
{
  LC_BEGIN_TESTS("lcobjcreator");
  test_1();
  LC_END_TESTS;
}
