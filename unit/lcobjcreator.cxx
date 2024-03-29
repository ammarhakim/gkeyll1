/**
 * @file	lccreator.cxx
 *
 * @brief	Unit tests for Lucee::ObjCreator and related classes
 */

// lib includes
#include <LcBasicObj.h>
#include <LcObjCreator.h>
#include <LcObjRegistry.h>
#include <LcTest.h>

class Base : public Lucee::BasicObj
{
  public:
    static const char *id;

    virtual ~Base() {}

    virtual std::string getMessage() = 0;

    virtual void readInput(Lucee::LuaTable& tbl)
    {
    }
};
const char *Base::id = "Base";

class Derived : public Base 
{
  public:
    static const char *id;
    std::string getMessage()
    {
      return "I am Derived";
    }
};
const char *Derived::id = "derived";

class Derived2 : public Base 
{
  public:
    static const char *id;
    std::string getMessage()
    {
      return "I am Derived2";
    }
};
const char *Derived2::id = "derived2";

Lucee::ObjRegistry<Base, Derived> derived;
Lucee::ObjRegistry<Base, Derived2> derived2;

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
main(int argc, char **argv)
{
  LC_BEGIN_TESTS("lcobjcreator");
  test_1();
  LC_END_TESTS;
}
