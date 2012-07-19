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
#include <valarray>

void
test_1()
{
  blitz::Array<double, 1> x(10), y(10), z(10);
  x = 1.0;
  y = 2.0;
  z = x+y;
}

void
test_2()
{
  blitz::Array<int,2> A(4,4);
  A = blitz::tensor::i * 10 + blitz::tensor::j; 
  std::cout << "A = " << A << std::endl;
  
  blitz::Array<int,2>::iterator iter = A.begin(), end = A.end();

  while (iter != end)
  {
    std::cout << iter.position() << '\t' << (*iter) << std::endl;
    ++iter;
  }
}

void
test_3()
{
  blitz::Array<blitz::TinyVector<float, 3>, 2> A(5,5);
  A = 0.0; // clear out contents
  blitz::Array<float, 2> B = A.extractComponent(float(), 1, 3);
  B = 10;
  std::cout << A << std::endl;
}

void
test_4()
{
  blitz::Array<blitz::TinyVector<float, 3>, 2> A(5,5);

  blitz::Array<float, 2> B0 = A.extractComponent(float(), 0, 1);
  blitz::Array<float, 2> B1 = A.extractComponent(float(), 1, 2);
  blitz::Array<float, 2> B2 = A.extractComponent(float(), 2, 3);
  B0 = 0;
  B1 = 1;
  B2 = 2;

  std::cout << A << std::endl;
}

int
main(int argc, char *argv[])
{
  test_1();
  test_2();
  test_3();
  test_4();
  return 0;
}
