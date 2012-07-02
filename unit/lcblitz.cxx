/**
 * @file	lcblitz.cxx
 *
 * @brief	Unit tests for using BLITZ in Lucee
 */

// blitz inlcudes
#include <blitz/array.h>

// std includes
#include <cstdlib>
#include <iostream>
#include <vector>

void
test_1()
{
  blitz::Array<double, 1> x(10), y(10), z(10);
  x = 1.0;
  y = 2.0;
  z = x+y;
}

int
main(int argc, char *argv[])
{
  test_1();
  return 0;
}
