/**
 * @file	lcfield.cxx
 *
 * @brief	Unit tests for Lucee::Field class
 */

// std includes
#include <iostream>
#include <cstring>

// lucee includes
#include <LcBase64.h>
#include <LcTest.h>

void
test_1()
{
  base64::encoder E;
  //E.encode(std::cin, std::cout);
  
}

int
main(void)
{
  LC_BEGIN_TESTS("lcbase64");
  test_1();
  LC_END_TESTS;
}


