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
    if (elmOrder == 2)
      this->setNumNodes(4);
    else if (elmOrder == 3)
      this->setNumNodes(8);
    else
    {
      Lucee::Except lce("SerendipityElement2D: Order must be 2 or 3.");
      lce << " Provided " << elmOrder << " instead";
      throw lce;
    }

// initialize matrices
    if (elmOrder == 2)
      setupOrder2();
    else if (elmOrder == 3)
      setupOrder3();
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
  SerendipityElement2D::setupOrder2()
  {
    unsigned shape[2] = {4,4};
    int start[2] = {1,1};

// mass matrix
    refNjNk = Lucee::Matrix<double>(shape, start);
    refNjNk(1,1) = 4.0/9.0;
    refNjNk(1,2) = 2.0/9.0;
    refNjNk(1,3) = 4.0/9.0;
    refNjNk(1,4) = 1.0/9.0;
    refNjNk(2,1) = 2.0/9.0;
    refNjNk(2,2) = 4.0/9.0;
    refNjNk(2,3) = 2.0/9.0;
    refNjNk(2,4) = 2.0/9.0;
    refNjNk(3,1) = 4.0/9.0;
    refNjNk(3,2) = 2.0/9.0;
    refNjNk(3,3) = 4.0/9.0;
    refNjNk(3,4) = 1.0/9.0;
    refNjNk(4,1) = 1.0/9.0;
    refNjNk(4,2) = 2.0/9.0;
    refNjNk(4,3) = 1.0/9.0;
    refNjNk(4,4) = 4.0/9.0;

// stiffness matrix
    refDNjDNk = Lucee::Matrix<double>(shape, start);
    refDNjDNk(1,1) = 2.0/3.0;
    refDNjDNk(1,2) = (-1.0)/6.0;
    refDNjDNk(1,3) = 2.0/3.0;
    refDNjDNk(1,4) = (-1.0)/3.0;
    refDNjDNk(2,1) = (-1.0)/6.0;
    refDNjDNk(2,2) = 2.0/3.0;
    refDNjDNk(2,3) = (-1.0)/6.0;
    refDNjDNk(2,4) = (-1.0)/6.0;
    refDNjDNk(3,1) = 2.0/3.0;
    refDNjDNk(3,2) = (-1.0)/6.0;
    refDNjDNk(3,3) = 2.0/3.0;
    refDNjDNk(3,4) = (-1.0)/3.0;
    refDNjDNk(4,1) = (-1.0)/3.0;
    refDNjDNk(4,2) = (-1.0)/6.0;
    refDNjDNk(4,3) = (-1.0)/3.0;
    refDNjDNk(4,4) = 2.0/3.0;
  }
}
