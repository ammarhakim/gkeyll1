/**
 * @file	lcfieldptr.cxx
 *
 * @brief	Unit tests for FieldPtr and ConstFieldPtr classes
 */

// lucee includes
#include <LcConstFieldPtr.h>
#include <LcFieldPtr.h>
#include <LcTest.h>

// std includes
#include <vector>

void
test_1(const std::vector<double>& x, Lucee::FieldPtr<double> fp)
{
  LC_ASSERT("Testing if size are the same", 
    x.size() == fp.getNumComponents());
  for (unsigned i=0; i<x.size(); ++i)
    LC_ASSERT("Testing if values are the same",
      x[i] == fp[i]);
}

void
setFp(Lucee::FieldPtr<double> fp)
{
  for (unsigned i=0; i<fp.getNumComponents(); ++i)
    fp[i] = (4.0+i)*0.5;
}

void
test_2(const std::vector<double>& x, const Lucee::ConstFieldPtr<double>& fp)
{
  LC_ASSERT("Testing if size are the same", 
    x.size() == fp.getNumComponents());
  for (unsigned i=0; i<x.size(); ++i)
    LC_ASSERT("Testing if values are the same",
      x[i] == fp[i]);
}

void
test_3()
{
  Lucee::FieldPtr<double> fp(10);
  fp = 3.5;
  for (unsigned i=0; i<fp.getNumComponents(); ++i)
    LC_ASSERT("Testing for assignment", fp[i] == 3.5);
  for (unsigned i=0; i<fp.getNumComponents(); ++i)
    fp[i] = (i+0.5)*10;
  for (unsigned i=0; i<fp.getNumComponents(); ++i)
    LC_ASSERT("Testing for assignment", fp[i] == (i+0.5)*10);
}

int
main(void) 
{
  LC_BEGIN_TESTS("lcfieldptr");
  std::vector<double> x(10);
  for (unsigned i=0; i<x.size(); ++i)
    x[i] = (i+0.5)*3.0;

  test_1(x, x);
  test_2(x, x);

  Lucee::FieldPtr<double> fp(x);
  test_1(x, fp);
  test_2(x, fp);

  Lucee::ConstFieldPtr<double> cfp(x);
  test_2(x, cfp);

  setFp(fp);
  for (unsigned i=0; i<x.size(); ++i)
    LC_ASSERT("Testing for set values", x[i] == (4.0+i)*0.5);

  std::vector<double> xx(10);
  setFp(xx);
  for (unsigned i=0; i<xx.size(); ++i)
    LC_ASSERT("Testing for set values", xx[i] == (4.0+i)*0.5);

  test_3();

  LC_END_TESTS;
}
