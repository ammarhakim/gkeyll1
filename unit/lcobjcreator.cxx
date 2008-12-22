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

// lib includes
#include <lcobjregister.h>
#include <lcobjcreator.h>
#include <lctest.h>

// etc includes
#include <loki/Factory.h>
#include <loki/Singleton.h>

template <typename B>
class FactoryBase
{
  public:
    typedef Loki::SingletonHolder<Loki::Factory<B, std::string> > ObjFactory;
};

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

Lucee::ObjRegister<Derived, Base> derived("derived");
Lucee::ObjRegister<Derived2, Base> derived2("derived2");

void
test_a()
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
  test_a();
  LC_END_TESTS;
}

