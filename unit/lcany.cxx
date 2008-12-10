/**
 * @file	lcany.cxx
 *
 * @brief	Unit tests for Lucee::Any class
 *
 * @version	$Id$ *
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#include <lcany.h>
#include <lctest.h>

#include <typeinfo>
#include <iostream>
#include <string>
#include <string.h>
#include <vector>

using namespace std;

void
test_any_a()
{ // test storing different value types
    
  Lucee::Any any; // create empty object

  LC_ASSERT("Testing if object is empty", any.isEmpty()==true);
  // insert int
  any = Lucee::Any(1);
  LC_ASSERT("Testing if object is empty", any.isEmpty()==false);    
  LC_ASSERT("Testing retrived value", Lucee::any_cast<int>(any)==1);

  // insert double
  any = Lucee::Any(3.1415);
  LC_ASSERT("Testing retrived value", Lucee::any_cast<double>(any)==3.1415);    

  // insert string
  any = Lucee::Any(string("Hello World"));
  LC_ASSERT("Testing retrived value", Lucee::any_cast<string>(any)==string("Hello World"));

  // insert char *
  any = Lucee::Any(strdup("Hello World"));
  LC_ASSERT("Testing retrived value", strcmp(Lucee::any_cast<char*>(any),"Hello World")==0);

  // insert double *
  double *x = new double[3];
  x[0] = 1; x[1] = 2; x[2] = 3;
  any = Lucee::Any(x);
  double *xx = Lucee::any_cast<double*>(any);
  for(unsigned i=0; i<3; ++i)
    LC_ASSERT("Testing retrived value", xx[i]==x[i]);

  // insert int and try removing it as double
  any = Lucee::Any(1);
  LC_RAISES("Removing double from int any", Lucee::any_cast<double>(any), bad_cast);

  // test copy ctor
  any = Lucee::Any(1);
  Lucee::Any any1;
  any1 = any;
  LC_ASSERT("Testing copied value", Lucee::any_cast<int>(any1)==1);

  // test swap
  any = Lucee::Any(1);
  any1 = Lucee::Any(2);
  any.swap(any1);
  LC_ASSERT("Testing swaped value", Lucee::any_cast<int>(any1)==1);
  LC_ASSERT("Testing swaped value", Lucee::any_cast<int>(any)==2);

  // test ass operator to a different type
  any1 = 2.0;
  LC_ASSERT("Testing swaped value", Lucee::any_cast<double>(any1)==2.0);

}

void
test_any_b()
{
  vector<Lucee::Any> v_wxa;

  v_wxa.push_back(22);
  v_wxa.push_back(string("Hello World"));
  v_wxa.push_back(3.1415);

  LC_ASSERT("Testing value 0", Lucee::any_cast<int>(v_wxa[0]) == 22);
  LC_ASSERT("Testing value 1", Lucee::any_cast<string>(v_wxa[1]) == "Hello World");
  LC_ASSERT("Testing value 2", Lucee::any_cast<double>(v_wxa[2]) == 3.1415);
}

void
test_any_c() 
{
  Lucee::Any iv(10);
  const void *ivp = iv.to_void_ptr();
  LC_ASSERT("Testing of void * of int Lucee::Any worked",
    *((int *)ivp) == 10);

  // insert double *
  double *x = new double[3];
  x[0] = 1; x[1] = 2; x[2] = 3;
  Lucee::Any dpv(x);
  const void* dpvp = dpv.to_void_ptr();
  const double *xx = *((const double**) dpvp);
  for(unsigned i=0; i<3; ++i)
    LC_ASSERT("Testing retrived value", xx[i]==x[i]);

}

int
main(void)
{
  LC_BEGIN_TESTS("lcany");
  test_any_a();
  test_any_b();
  test_any_c();
  LC_END_TESTS;
}
