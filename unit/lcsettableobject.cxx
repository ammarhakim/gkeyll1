/**
 * @file	lcsettableobject.cxx
 *
 * @brief	Unit tests for Lucee::SettableObject class
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// lucee includes
#include <LcTest.h>
#include <LcSettableObject.h>
#include <LcExcept.h>

class Grid : public Lucee::SettableObject
{
  public:
    Grid()
      : Lucee::SettableObject("grid")
    {
      addData("numCells", &ncell, "Number of cells in grid");
      addData("length", &xlen, "Length of domain");
      addData("vec", &vec, "Vector in domain");
      addData("name", &name, "Name of grid");
    }

    bool init()
    {
      return true;
    }

    int ncell;
    double xlen;
    std::string name;
    std::vector<int> vec;
};

void
test_1()
{
  Grid grid;

  LC_RAISES("Testing to see if incorrect data can be set", 
    grid.setData("foo", 3), Lucee::Except);

  grid.setData("numCells", 10);
  grid.setData("length", 1.0);
  grid.setData("name", std::string("cart1D"));

  std::vector<int> myVec(3);
  myVec[0] = 1; myVec[1] = 2; myVec[2] = 3;
  grid.setData("vec", myVec);

  LC_ASSERT("Checking if setData worked", grid.ncell==10);
  LC_ASSERT("Checking if setData worked", grid.xlen==1.0);
  LC_ASSERT("Checking if setData worked", grid.name=="cart1D");

  LC_ASSERT("Checking if setData worked", grid.vec[0] == 1);
  LC_ASSERT("Checking if setData worked", grid.vec[1] == 2);
  LC_ASSERT("Checking if setData worked", grid.vec[2] == 3);
}

int
main(void)
{
  LC_BEGIN_TESTS("lcsettableobject");
  test_1();
  LC_END_TESTS;
}
