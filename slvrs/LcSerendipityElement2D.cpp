/**
 * @file	LcSerendipityElement2D.cpp
 *
 * @brief       Reference finite element with serendipity basis
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSerendipityElement2D.h>

namespace Lucee
{
  const char *SerendipityElement2D::id = "Serendipity2D";

  SerendipityElement2D::SerendipityElement2D(unsigned order)
    : Lucee::NodalFiniteElementIfc(4), elmOrder(order), 
      refNjNk((unsigned)0,(unsigned)0), refDNjDNk((unsigned)0,(unsigned)0)
  {
    if (elmOrder == 1)
      this->setNumNodes(4);
    else if (elmOrder == 2)
      this->setNumNodes(8);
    else
    {
      Lucee::Except lce("SerendipityElement2D: Order must be 1 or 2.");
      lce << " Provided " << elmOrder << " instead";
      throw lce;
    }

// initialize matrices
    if (elmOrder == 1)
      setupOrder1();
    else if (elmOrder == 2)
      setupOrder2();
  }

  void
  SerendipityElement2D::getMassMatrix(Lucee::Matrix<double> NjNk) const
  {
    NjNk = refNjNk;
  }

  void
  SerendipityElement2D::getStiffnessMatrix(Lucee::Matrix<double> DNjDNk) const
  {
    DNjDNk = refDNjDNk;
  }

  void
  SerendipityElement2D::setupOrder1()
  {
    unsigned shape[2] = {4,4};
    int start[2] = {1,1};

// mass matrix (automatically generated. See scripts/serendipity-2D.mac)
    refNjNk = Lucee::Matrix<double>(shape, start);
    refNjNk(1,1) = 4.0/9.0;
    refNjNk(1,2) = 2.0/9.0;
    refNjNk(1,3) = 1.0/9.0;
    refNjNk(1,4) = 2.0/9.0;
    refNjNk(2,1) = 2.0/9.0;
    refNjNk(2,2) = 4.0/9.0;
    refNjNk(2,3) = 2.0/9.0;
    refNjNk(2,4) = 1.0/9.0;
    refNjNk(3,1) = 1.0/9.0;
    refNjNk(3,2) = 2.0/9.0;
    refNjNk(3,3) = 4.0/9.0;
    refNjNk(3,4) = 2.0/9.0;
    refNjNk(4,1) = 2.0/9.0;
    refNjNk(4,2) = 1.0/9.0;
    refNjNk(4,3) = 2.0/9.0;
    refNjNk(4,4) = 4.0/9.0;

// stiffness matrix (automatically generated. See scripts/serendipity-2D.mac)
    refDNjDNk = Lucee::Matrix<double>(shape, start);
    refDNjDNk(1,1) = 2.0/3.0;
    refDNjDNk(1,2) = (-1.0)/6.0;
    refDNjDNk(1,3) = (-1.0)/3.0;
    refDNjDNk(1,4) = (-1.0)/6.0;
    refDNjDNk(2,1) = (-1.0)/6.0;
    refDNjDNk(2,2) = 2.0/3.0;
    refDNjDNk(2,3) = (-1.0)/6.0;
    refDNjDNk(2,4) = (-1.0)/3.0;
    refDNjDNk(3,1) = (-1.0)/3.0;
    refDNjDNk(3,2) = (-1.0)/6.0;
    refDNjDNk(3,3) = 2.0/3.0;
    refDNjDNk(3,4) = (-1.0)/6.0;
    refDNjDNk(4,1) = (-1.0)/6.0;
    refDNjDNk(4,2) = (-1.0)/3.0;
    refDNjDNk(4,3) = (-1.0)/6.0;
    refDNjDNk(4,4) = 2.0/3.0;
  }

  void
  SerendipityElement2D::setupOrder2()
  {
    unsigned shape[2] = {8,8};
    int start[2] = {1,1};

// mass matrix (automatically generated. See scripts/serendipity-2D.mac)
    refNjNk = Lucee::Matrix<double>(shape, start);
    refNjNk(1,1) = 2.0/15.0;
    refNjNk(1,2) = 2.0/45.0;
    refNjNk(1,3) = 1.0/15.0;
    refNjNk(1,4) = 2.0/45.0;
    refNjNk(1,5) = (-2.0)/15.0;
    refNjNk(1,6) = (-8.0)/45.0;
    refNjNk(1,7) = (-8.0)/45.0;
    refNjNk(1,8) = (-2.0)/15.0;
    refNjNk(2,1) = 2.0/45.0;
    refNjNk(2,2) = 2.0/15.0;
    refNjNk(2,3) = 2.0/45.0;
    refNjNk(2,4) = 1.0/15.0;
    refNjNk(2,5) = (-2.0)/15.0;
    refNjNk(2,6) = (-2.0)/15.0;
    refNjNk(2,7) = (-8.0)/45.0;
    refNjNk(2,8) = (-8.0)/45.0;
    refNjNk(3,1) = 1.0/15.0;
    refNjNk(3,2) = 2.0/45.0;
    refNjNk(3,3) = 2.0/15.0;
    refNjNk(3,4) = 2.0/45.0;
    refNjNk(3,5) = (-8.0)/45.0;
    refNjNk(3,6) = (-2.0)/15.0;
    refNjNk(3,7) = (-2.0)/15.0;
    refNjNk(3,8) = (-8.0)/45.0;
    refNjNk(4,1) = 2.0/45.0;
    refNjNk(4,2) = 1.0/15.0;
    refNjNk(4,3) = 2.0/45.0;
    refNjNk(4,4) = 2.0/15.0;
    refNjNk(4,5) = (-8.0)/45.0;
    refNjNk(4,6) = (-8.0)/45.0;
    refNjNk(4,7) = (-2.0)/15.0;
    refNjNk(4,8) = (-2.0)/15.0;
    refNjNk(5,1) = (-2.0)/15.0;
    refNjNk(5,2) = (-2.0)/15.0;
    refNjNk(5,3) = (-8.0)/45.0;
    refNjNk(5,4) = (-8.0)/45.0;
    refNjNk(5,5) = 32.0/45.0;
    refNjNk(5,6) = 4.0/9.0;
    refNjNk(5,7) = 16.0/45.0;
    refNjNk(5,8) = 4.0/9.0;
    refNjNk(6,1) = (-8.0)/45.0;
    refNjNk(6,2) = (-2.0)/15.0;
    refNjNk(6,3) = (-2.0)/15.0;
    refNjNk(6,4) = (-8.0)/45.0;
    refNjNk(6,5) = 4.0/9.0;
    refNjNk(6,6) = 32.0/45.0;
    refNjNk(6,7) = 4.0/9.0;
    refNjNk(6,8) = 16.0/45.0;
    refNjNk(7,1) = (-8.0)/45.0;
    refNjNk(7,2) = (-8.0)/45.0;
    refNjNk(7,3) = (-2.0)/15.0;
    refNjNk(7,4) = (-2.0)/15.0;
    refNjNk(7,5) = 16.0/45.0;
    refNjNk(7,6) = 4.0/9.0;
    refNjNk(7,7) = 32.0/45.0;
    refNjNk(7,8) = 4.0/9.0;
    refNjNk(8,1) = (-2.0)/15.0;
    refNjNk(8,2) = (-8.0)/45.0;
    refNjNk(8,3) = (-8.0)/45.0;
    refNjNk(8,4) = (-2.0)/15.0;
    refNjNk(8,5) = 4.0/9.0;
    refNjNk(8,6) = 16.0/45.0;
    refNjNk(8,7) = 4.0/9.0;
    refNjNk(8,8) = 32.0/45.0;

// stiffness matrix (automatically generated. See scripts/serendipity-2D.mac)
    refDNjDNk = Lucee::Matrix<double>(shape, start);
    refDNjDNk(1,1) = 52.0/45.0;
    refDNjDNk(1,2) = 1.0/2.0;
    refDNjDNk(1,3) = 23.0/45.0;
    refDNjDNk(1,4) = 1.0/2.0;
    refDNjDNk(1,5) = (-37.0)/45.0;
    refDNjDNk(1,6) = (-23.0)/45.0;
    refDNjDNk(1,7) = (-23.0)/45.0;
    refDNjDNk(1,8) = (-37.0)/45.0;
    refDNjDNk(2,1) = 1.0/2.0;
    refDNjDNk(2,2) = 52.0/45.0;
    refDNjDNk(2,3) = 1.0/2.0;
    refDNjDNk(2,4) = 23.0/45.0;
    refDNjDNk(2,5) = (-37.0)/45.0;
    refDNjDNk(2,6) = (-37.0)/45.0;
    refDNjDNk(2,7) = (-23.0)/45.0;
    refDNjDNk(2,8) = (-23.0)/45.0;
    refDNjDNk(3,1) = 23.0/45.0;
    refDNjDNk(3,2) = 1.0/2.0;
    refDNjDNk(3,3) = 52.0/45.0;
    refDNjDNk(3,4) = 1.0/2.0;
    refDNjDNk(3,5) = (-23.0)/45.0;
    refDNjDNk(3,6) = (-37.0)/45.0;
    refDNjDNk(3,7) = (-37.0)/45.0;
    refDNjDNk(3,8) = (-23.0)/45.0;
    refDNjDNk(4,1) = 1.0/2.0;
    refDNjDNk(4,2) = 23.0/45.0;
    refDNjDNk(4,3) = 1.0/2.0;
    refDNjDNk(4,4) = 52.0/45.0;
    refDNjDNk(4,5) = (-23.0)/45.0;
    refDNjDNk(4,6) = (-23.0)/45.0;
    refDNjDNk(4,7) = (-37.0)/45.0;
    refDNjDNk(4,8) = (-37.0)/45.0;
    refDNjDNk(5,1) = (-37.0)/45.0;
    refDNjDNk(5,2) = (-37.0)/45.0;
    refDNjDNk(5,3) = (-23.0)/45.0;
    refDNjDNk(5,4) = (-23.0)/45.0;
    refDNjDNk(5,5) = 104.0/45.0;
    refDNjDNk(5,6) = 0;
    refDNjDNk(5,7) = 16.0/45.0;
    refDNjDNk(5,8) = 0;
    refDNjDNk(6,1) = (-23.0)/45.0;
    refDNjDNk(6,2) = (-37.0)/45.0;
    refDNjDNk(6,3) = (-37.0)/45.0;
    refDNjDNk(6,4) = (-23.0)/45.0;
    refDNjDNk(6,5) = 0;
    refDNjDNk(6,6) = 104.0/45.0;
    refDNjDNk(6,7) = 0;
    refDNjDNk(6,8) = 16.0/45.0;
    refDNjDNk(7,1) = (-23.0)/45.0;
    refDNjDNk(7,2) = (-23.0)/45.0;
    refDNjDNk(7,3) = (-37.0)/45.0;
    refDNjDNk(7,4) = (-37.0)/45.0;
    refDNjDNk(7,5) = 16.0/45.0;
    refDNjDNk(7,6) = 0;
    refDNjDNk(7,7) = 104.0/45.0;
    refDNjDNk(7,8) = 0;
    refDNjDNk(8,1) = (-37.0)/45.0;
    refDNjDNk(8,2) = (-23.0)/45.0;
    refDNjDNk(8,3) = (-23.0)/45.0;
    refDNjDNk(8,4) = (-37.0)/45.0;
    refDNjDNk(8,5) = 0;
    refDNjDNk(8,6) = 16.0/45.0;
    refDNjDNk(8,7) = 0;
    refDNjDNk(8,8) = 104.0/45.0;
  }
}
