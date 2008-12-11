/**
 * @file	lckeyval.cxx
 *
 * @brief	Unit tests for Lucee::KeyVal class
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#include <lckeyval.h>
#include <lctest.h>
#include <lcexcept.h>

#include <iostream>
#include <string>
#include <stdexcept>

using namespace std;

void
test_a_kv(const Lucee::KeyVal& kv)
{
  LC_ASSERT("Testing retreving int", kv.get<int>("meqn")==18);
  LC_ASSERT("Testing retreving double", kv.get<double>("gas_gamma")==1.4);
  LC_ASSERT("Testing if algorithm exists", kv.has("algorithm") == true);
  LC_ASSERT("Testing retreving string", kv.get<string>("algorithm") == string("RKDG"));

  LC_ASSERT("Testing if humbug exists", kv.has("Humbug") == false);
  LC_RAISES("Testing if runtime is raised", kv.get<int>("Humbug"), Lucee::Except);
  LC_RAISES("Testing if bad_cast is raised", kv.get<double>("meqn"), bad_cast);


  std::vector<std::string>::const_iterator itr;
  std::vector<std::string> names;

  names = kv.getKeys<int>();
  LC_ASSERT("Testing name of first integer", names[0] == "meqn");
  LC_ASSERT("Testing name of first integer", names[1] == "int_a");
  LC_ASSERT("Testing name of second integer", names[2] == "int_b");
  LC_ASSERT("Testing name of third integer", names[3] == "int_c");

  names = kv.getKeys<double>();
  LC_ASSERT("Testing name of first double", names[0] == "gas_gamma");
  LC_ASSERT("Testing name of first double", names[1] == "double_a");
  LC_ASSERT("Testing name of second double", names[2] == "double_b");
  LC_ASSERT("Testing name of third double", names[3] == "double_c");

  std::vector<int> newValues 
    = kv.get<std::vector<int> >("intList");
  LC_ASSERT("Testing if size of returned array is correct",
    newValues.size() == 3);
  LC_ASSERT("Testing if elements of returned array are correct",
    newValues[0] == 1);
  LC_ASSERT("Testing if elements of returned array are correct",
    newValues[1] == 2);
  LC_ASSERT("Testing if elements of returned array are correct",
    newValues[2] == 3);

  names = kv.getKeys<std::vector<int> >();
  LC_ASSERT("Testing if names are proper",
    names[0] == "intList");
}

void
test_lckv_a()
{
  Lucee::KeyVal kv;

  kv.add("meqn", 18);
  kv.add("gas_gamma", 1.4);
  kv.add("algorithm", string("RKDG"));

  LC_ASSERT("Testing retreving int", kv.get<int>("meqn")==18);
  LC_ASSERT("Testing retreving double", kv.get<double>("gas_gamma")==1.4);
  LC_ASSERT("Testing retreving string", kv.get<string>("algorithm") == string("RKDG"));

  LC_RAISES("Testing if runtime is raised", kv.get<int>("Humbug"), Lucee::Except);
  LC_RAISES("Testing if bad_cast is raised", kv.get<double>("meqn"), bad_cast);

  kv.add("int_a", 1);
  kv.add("int_b", 2);
  kv.add("int_c", 3);

  kv.add("double_a", 1.0);
  kv.add("double_b", 2.0);
  kv.add("double_c", 3.0);

  std::vector<std::string>::const_iterator itr;
  std::vector<std::string> names;

  names = kv.getKeys<int>();
  LC_ASSERT("Testing name of first integer", names[0] == "meqn");
  LC_ASSERT("Testing name of first integer", names[1] == "int_a");
  LC_ASSERT("Testing name of second integer", names[2] == "int_b");
  LC_ASSERT("Testing name of third integer", names[3] == "int_c");

  names = kv.getKeys<double>();
  LC_ASSERT("Testing name of first double", names[0] == "gas_gamma");
  LC_ASSERT("Testing name of first double", names[1] == "double_a");
  LC_ASSERT("Testing name of second double", names[2] == "double_b");
  LC_ASSERT("Testing name of third double", names[3] == "double_c");

  std::vector<int> values;
  values.push_back(1);
  values.push_back(2);
  values.push_back(3);

  kv.add("intList", values);

  std::vector<int> newValues 
    = kv.get<std::vector<int> >("intList");
  LC_ASSERT("Testing if size of returned array is correct",
    newValues.size() == 3);
  LC_ASSERT("Testing if elements of returned array are correct",
    newValues[0] == 1);
  LC_ASSERT("Testing if elements of returned array are correct",
    newValues[1] == 2);
  LC_ASSERT("Testing if elements of returned array are correct",
    newValues[2] == 3);

  names = kv.getKeys<std::vector<int> >();
  LC_ASSERT("Testing if names are proper",
    names[0] == "intList");

  // make a copy
  Lucee::KeyVal kv_a(kv);
  // test if it is properly made
  test_a_kv(kv_a);

  Lucee::KeyVal kv_b;
  kv_b = kv_a;
  test_a_kv(kv_b);
}

int
main(void)
{
  LC_BEGIN_TESTS("lckeyval");
  test_lckv_a();
  LC_END_TESTS;
}
