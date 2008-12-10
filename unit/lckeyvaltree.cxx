/**
 * @file	lckeyvaltree.cxx
 *
 * @brief	Unit tests for Lucee::KeyValTree class
 *
 * @version	$Id$ *
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */
#include <lckeyvaltree.h>
#include <lctest.h>
#include <lcexcept.h>

#include <iostream>
#include <vector>
#include <string>

using namespace std;

void
test_a_keyvaltree(const Lucee::KeyValTree& kvt)
{
  LC_ASSERT("Testing int", kvt.get<int>("meqn")==6);
  LC_ASSERT("Testing double", kvt.get<double>("gas_gamma")==1.4);
  LC_ASSERT("Testing string", kvt.get<string>("algo")==string("WAVE"));

  LC_RAISES("Testing if wrong data can be gotten", kvt.get<double>("gasgamma"), Lucee::Except);
  LC_RAISES("Testing if wrong type can be gotten", kvt.get<int>("gas_gamma"), bad_cast);
}

void
test_keyvaltree_a()
{
  Lucee::KeyValTree kvt("Fluid");

  kvt.add("gas_gamma", 1.4);
  kvt.add("meqn", 6);
  kvt.add("algo", string("WAVE"));

  test_a_keyvaltree(kvt);

  // test copy ctor
  Lucee::KeyValTree kvt_a(kvt);
  test_a_keyvaltree(kvt_a); 

  // test ass operator
  kvt_a = kvt;
  test_a_keyvaltree(kvt_a);
}

void
test_keyvaltree_b()
{
  // top level set
  Lucee::KeyValTree kvt("warpx");
  //kvt.add("description", "Two-fluid solver");

  // ion set
  Lucee::KeyValTree kvt_ion("ions");
  kvt_ion.add("gas_gamma", 1.4);
  kvt_ion.add("meqn", 5);
  // .. add ion set
  kvt.addSet(kvt_ion);

  // electron set
  Lucee::KeyValTree kvt_elc("electrons");
  kvt_elc.add("gas_gamma", 2.0);
  kvt_elc.add("meqn", 5);

  // properties in electron set
  Lucee::KeyValTree kvt_elc_prop("properties");
  kvt_elc_prop.add("mass", 1.0);
  kvt_elc.addSet(kvt_elc_prop);

  kvt.addSet(kvt_elc);

  LC_ASSERT("Testing if ions keyvaltree exists", kvt.hasSet("ions")==true);
  const Lucee::KeyValTree& kvtcr = kvt.getSet("ions");

  LC_ASSERT("Testing for meqn", kvtcr.get<int>("meqn")==5);
  LC_ASSERT("Testing for gas_gamma", kvtcr.get<double>("gas_gamma")==1.4);

  const Lucee::KeyValTree& kvtcra = kvt.getSet("electrons");

  LC_ASSERT("Testing for meqn", kvtcra.get<int>("meqn")==5);
  LC_ASSERT("Testing for gas_gamma", kvtcra.get<double>("gas_gamma")==2.0);

  const Lucee::KeyValTree& kvtcrab = kvtcra.getSet("properties");
  LC_ASSERT("Testing for mass", kvtcrab.get<double>("mass")==1.0);

  LC_RAISES("Testing if bad set can be gotten", kvt.getSet("Flu"), Lucee::Except);
}

void
test_keyvaltree_d()
{
  Lucee::KeyValTree kvt("warpx");

  LC_ASSERT("Testing if name is correct",
    kvt.getName() == "warpx");

  Lucee::KeyValTree kvta("warpx/a");
  kvta.add<std::string>("Type", "type_a");
  Lucee::KeyValTree kvtb("warpx/b");
  kvtb.add<std::string>("Type", "type_a");
  Lucee::KeyValTree kvtc("warpx/c");
  kvtc.add<std::string>("Type", "type_a");

  kvt.addSet(kvta);
  kvt.addSet(kvtb);
  kvt.addSet(kvtc);

  std::vector<std::string> names = 
    kvt.getNamesOfType("type_a");
  LC_ASSERT("Testing names of type list is of correct size",
    names.size() == 3);
}

void
test_keyvaltree_e()
{
  Lucee::KeyValTree kvtv("warpx");

  kvtv.add("int_a", 1);
  kvtv.add("int_b", 2);
  kvtv.add("int_c", 3);

  kvtv.add("double_a", 1.0);
  kvtv.add("double_b", 2.0);
  kvtv.add("double_c", 3.0);

  std::vector<std::string>::const_iterator itr;
  std::vector<std::string> names;

  names = kvtv.getKeys<int>();
  LC_ASSERT("Testing name of first integer", names[0] == "int_a");
  LC_ASSERT("Testing name of second integer", names[1] == "int_b");
  LC_ASSERT("Testing name of third integer", names[2] == "int_c");

  names = kvtv.getKeys<double>();
  LC_ASSERT("Testing name of first double", names[0] == "double_a");
  LC_ASSERT("Testing name of second double", names[1] == "double_b");
  LC_ASSERT("Testing name of third double", names[2] == "double_c");

  Lucee::KeyValTree kvtvc = kvtv;

  names = kvtvc.getKeys<int>();
  LC_ASSERT("Testing name of first integer", names[0] == "int_a");
  LC_ASSERT("Testing name of second integer", names[1] == "int_b");
  LC_ASSERT("Testing name of third integer", names[2] == "int_c");

  names = kvtvc.getKeys<double>();
  LC_ASSERT("Testing name of first double", names[0] == "double_a");
  LC_ASSERT("Testing name of second double", names[1] == "double_b");
  LC_ASSERT("Testing name of third double", names[2] == "double_c");

  Lucee::KeyValTree kvtva;
  kvtva = kvtvc;

  names = kvtva.getKeys<int>();
  LC_ASSERT("Testing name of first integer", names[0] == "int_a");
  LC_ASSERT("Testing name of second integer", names[1] == "int_b");
  LC_ASSERT("Testing name of third integer", names[2] == "int_c");

  names = kvtva.getKeys<double>();
  LC_ASSERT("Testing name of first double", names[0] == "double_a");
  LC_ASSERT("Testing name of second double", names[1] == "double_b");
  LC_ASSERT("Testing name of third double", names[2] == "double_c");
}

int
main(void)
{
  LC_BEGIN_TESTS("lckeyvaltree");
  test_keyvaltree_a();
  test_keyvaltree_b();
  test_keyvaltree_d();
  test_keyvaltree_e();
  LC_END_TESTS;
}
