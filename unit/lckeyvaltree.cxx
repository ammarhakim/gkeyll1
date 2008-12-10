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
  Lucee::KeyValTree wx("warpx");
  //wx.add("description", "Two-fluid solver");

  // ion set
  Lucee::KeyValTree wx_ion("ions");
  wx_ion.add("gas_gamma", 1.4);
  wx_ion.add("meqn", 5);
  // .. add ion set
  wx.addSet(wx_ion);

  // electron set
  Lucee::KeyValTree wx_elc("electrons");
  wx_elc.add("gas_gamma", 2.0);
  wx_elc.add("meqn", 5);

  // properties in electron set
  Lucee::KeyValTree wx_elc_prop("properties");
  wx_elc_prop.add("mass", 1.0);
  wx_elc.addSet(wx_elc_prop);

  wx.addSet(wx_elc);

  LC_ASSERT("Testing if ions keyvaltree exists", wx.hasSet("ions")==true);
  const Lucee::KeyValTree& wxcr = wx.getSet("ions");

  LC_ASSERT("Testing for meqn", wxcr.get<int>("meqn")==5);
  LC_ASSERT("Testing for gas_gamma", wxcr.get<double>("gas_gamma")==1.4);

  const Lucee::KeyValTree& wxcra = wx.getSet("electrons");

  LC_ASSERT("Testing for meqn", wxcra.get<int>("meqn")==5);
  LC_ASSERT("Testing for gas_gamma", wxcra.get<double>("gas_gamma")==2.0);

  const Lucee::KeyValTree& wxcrab = wxcra.getSet("properties");
  LC_ASSERT("Testing for mass", wxcrab.get<double>("mass")==1.0);

  LC_RAISES("Testing if bad set can be gotten", wx.getSet("Flu"), Lucee::Except);
}

void
test_keyvaltree_d()
{
  Lucee::KeyValTree wx("warpx");

  Lucee::KeyValTree wxa("warpx/a");
  wxa.add<std::string>("Type", "type_a");
  Lucee::KeyValTree wxb("warpx/b");
  wxb.add<std::string>("Type", "type_a");
  Lucee::KeyValTree wxc("warpx/c");
  wxc.add<std::string>("Type", "type_a");

  wx.addSet(wxa);
  wx.addSet(wxb);
  wx.addSet(wxc);

  std::vector<std::string> names = 
    wx.getNamesOfType("type_a");
  LC_ASSERT("Testing names of type list is of correct size",
    names.size() == 3);
}

void
test_keyvaltree_e()
{
  Lucee::KeyValTree wxv("warpx");

  wxv.add("int_a", 1);
  wxv.add("int_b", 2);
  wxv.add("int_c", 3);

  wxv.add("double_a", 1.0);
  wxv.add("double_b", 2.0);
  wxv.add("double_c", 3.0);

  std::vector<std::string>::const_iterator itr;
  std::vector<std::string> names;

  names = wxv.getKeys<int>();
  LC_ASSERT("Testing name of first integer", names[0] == "int_a");
  LC_ASSERT("Testing name of second integer", names[1] == "int_b");
  LC_ASSERT("Testing name of third integer", names[2] == "int_c");

  names = wxv.getKeys<double>();
  LC_ASSERT("Testing name of first double", names[0] == "double_a");
  LC_ASSERT("Testing name of second double", names[1] == "double_b");
  LC_ASSERT("Testing name of third double", names[2] == "double_c");

  Lucee::KeyValTree wxvc = wxv;

  names = wxvc.getKeys<int>();
  LC_ASSERT("Testing name of first integer", names[0] == "int_a");
  LC_ASSERT("Testing name of second integer", names[1] == "int_b");
  LC_ASSERT("Testing name of third integer", names[2] == "int_c");

  names = wxvc.getKeys<double>();
  LC_ASSERT("Testing name of first double", names[0] == "double_a");
  LC_ASSERT("Testing name of second double", names[1] == "double_b");
  LC_ASSERT("Testing name of third double", names[2] == "double_c");

  Lucee::KeyValTree wxva;
  wxva = wxvc;

  names = wxva.getKeys<int>();
  LC_ASSERT("Testing name of first integer", names[0] == "int_a");
  LC_ASSERT("Testing name of second integer", names[1] == "int_b");
  LC_ASSERT("Testing name of third integer", names[2] == "int_c");

  names = wxva.getKeys<double>();
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
