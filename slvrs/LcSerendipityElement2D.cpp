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
#include <LcStructuredGridBase.h>

namespace Lucee
{
  const char *SerendipityElement2D::id = "Serendipity2D";

  SerendipityElement2D::SerendipityElement2D()
    : Lucee::NodalFiniteElementIfc(4), polyOrder(2), 
      refNjNk(polyOrder,polyOrder), refDNjDNk(polyOrder,polyOrder)
  {
  }

  void
  SerendipityElement2D::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::NodalFiniteElementIfc::readInput(tbl);

// read in polynomial order
    polyOrder = tbl.getNumber("polyOrder");

    if (polyOrder == 1)
      this->setNumNodes(4);
    else if (polyOrder == 2)
      this->setNumNodes(8);
    else
    {
      Lucee::Except lce("SerendipityElement2D: Order must be 1 or 2.");
      lce << " Provided " << polyOrder << " instead";
      throw lce;
    }

// initialize matrices
    if (polyOrder == 1)
      setupOrder1();
    else if (polyOrder == 2)
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

// get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();

// get grid spacing (this is assumed to be uniform for now)
    double dx = grid.getDx(0), dy = grid.getDx(1);
    double dx2 = dx*dx, dy2 = dy*dy;

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

// scale to bring this into physical space
    refNjNk *= 0.5*dx*0.5*dy;

// stiffness matrix (automatically generated. See scripts/serendipity-2D.mac)
    refDNjDNk = Lucee::Matrix<double>(shape, start);
    refDNjDNk(1,1) = (7*dy2+4*dx2)/(dx2*dy2)/6.0+(dy2+4*dx2)/(dx2*dy2)/6.0;
    refDNjDNk(1,2) = -(8*dy2-4*dx2)/(dx2*dy2)/6.0;
    refDNjDNk(1,3) = (-4*dy2-4*dx2)/(dx2*dy2)/6.0;
    refDNjDNk(1,4) = 2.0*(dy2-2*dx2)/(3.0*dx2*dy2);
    refDNjDNk(2,1) = -(8*dy2-4*dx2)/(dx2*dy2)/6.0;
    refDNjDNk(2,2) = (7*dy2+4*dx2)/(dx2*dy2)/6.0+(dy2+4*dx2)/(dx2*dy2)/6.0;
    refDNjDNk(2,3) = 2.0*(dy2-2*dx2)/(3.0*dx2*dy2);
    refDNjDNk(2,4) = (-4*dy2-4*dx2)/(dx2*dy2)/6.0;
    refDNjDNk(3,1) = (-4*dy2-4*dx2)/(dx2*dy2)/6.0;
    refDNjDNk(3,2) = 2.0*(dy2-2*dx2)/(3.0*dx2*dy2);
    refDNjDNk(3,3) = (7*dy2+4*dx2)/(dx2*dy2)/6.0+(dy2+4*dx2)/(dx2*dy2)/6.0;
    refDNjDNk(3,4) = -(8*dy2-4*dx2)/(dx2*dy2)/6.0;
    refDNjDNk(4,1) = 2.0*(dy2-2*dx2)/(3.0*dx2*dy2);
    refDNjDNk(4,2) = (-4*dy2-4*dx2)/(dx2*dy2)/6.0;
    refDNjDNk(4,3) = -(8*dy2-4*dx2)/(dx2*dy2)/6.0;
    refDNjDNk(4,4) = (7*dy2+4*dx2)/(dx2*dy2)/6.0+(dy2+4*dx2)/(dx2*dy2)/6.0;

// scale to bring this into physical space
    refDNjDNk *= 0.5*dx*0.5*dy;
  }

  void
  SerendipityElement2D::setupOrder2()
  {
    unsigned shape[2] = {8,8};
    int start[2] = {1,1};

// get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();

// get grid spacing (this is assumed to be uniform for now)
    double dx = grid.getDx(0), dy = grid.getDx(1);
    double dx2 = dx*dx, dy2 = dy*dy;

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

// scale to bring this into physical space
    refNjNk *= 0.5*dx*0.5*dy;

// stiffness matrix (automatically generated. See scripts/serendipity-2D.mac)
    refDNjDNk = Lucee::Matrix<double>(shape, start);
    refDNjDNk(1,1) = (373*dy2+328*dx2)/(dx2*dy2)/180.0+(43*dy2+88*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(1,2) = -(-(187*dy2+68*dx2)/6.0-(37*dy2+68*dx2)/6.0)/(dx2*dy2)/30.0;
    refDNjDNk(1,3) = 2.0*(23*dy2+23*dx2)/(45.0*dx2*dy2);
    refDNjDNk(1,4) = 2.0*(17*dy2+28*dx2)/(45.0*dx2*dy2);
    refDNjDNk(1,5) = -(140*dy2-36*dx2)/(dx2*dy2)/45.0-(20*dy2+24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(1,6) = (21*dy2-160*dx2)/(dx2*dy2)/180.0-(69*dy2+160*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(1,7) = -(40*dy2+36*dx2)/(dx2*dy2)/45.0-(40*dy2-24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(1,8) = (69*dy2-440*dx2)/(dx2*dy2)/180.0-(21*dy2+200*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(2,1) = -(-(187*dy2+68*dx2)/6.0-(37*dy2+68*dx2)/6.0)/(dx2*dy2)/30.0;
    refDNjDNk(2,2) = (373*dy2+328*dx2)/(dx2*dy2)/180.0+(43*dy2+88*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(2,3) = 2.0*(17*dy2+28*dx2)/(45.0*dx2*dy2);
    refDNjDNk(2,4) = 2.0*(23*dy2+23*dx2)/(45.0*dx2*dy2);
    refDNjDNk(2,5) = -(140*dy2-36*dx2)/(dx2*dy2)/45.0-(20*dy2+24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(2,6) = (69*dy2-440*dx2)/(dx2*dy2)/180.0-(21*dy2+200*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(2,7) = -(40*dy2+36*dx2)/(dx2*dy2)/45.0-(40*dy2-24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(2,8) = (21*dy2-160*dx2)/(dx2*dy2)/180.0-(69*dy2+160*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(3,1) = 2.0*(23*dy2+23*dx2)/(45.0*dx2*dy2);
    refDNjDNk(3,2) = 2.0*(17*dy2+28*dx2)/(45.0*dx2*dy2);
    refDNjDNk(3,3) = (373*dy2+328*dx2)/(dx2*dy2)/180.0+(43*dy2+88*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(3,4) = -(-(187*dy2+68*dx2)/6.0-(37*dy2+68*dx2)/6.0)/(dx2*dy2)/30.0;
    refDNjDNk(3,5) = -(40*dy2+36*dx2)/(dx2*dy2)/45.0-(40*dy2-24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(3,6) = (69*dy2-440*dx2)/(dx2*dy2)/180.0-(21*dy2+200*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(3,7) = -(140*dy2-36*dx2)/(dx2*dy2)/45.0-(20*dy2+24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(3,8) = (21*dy2-160*dx2)/(dx2*dy2)/180.0-(69*dy2+160*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(4,1) = 2.0*(17*dy2+28*dx2)/(45.0*dx2*dy2);
    refDNjDNk(4,2) = 2.0*(23*dy2+23*dx2)/(45.0*dx2*dy2);
    refDNjDNk(4,3) = -(-(187*dy2+68*dx2)/6.0-(37*dy2+68*dx2)/6.0)/(dx2*dy2)/30.0;
    refDNjDNk(4,4) = (373*dy2+328*dx2)/(dx2*dy2)/180.0+(43*dy2+88*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(4,5) = -(40*dy2+36*dx2)/(dx2*dy2)/45.0-(40*dy2-24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(4,6) = (21*dy2-160*dx2)/(dx2*dy2)/180.0-(69*dy2+160*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(4,7) = -(140*dy2-36*dx2)/(dx2*dy2)/45.0-(20*dy2+24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(4,8) = (69*dy2-440*dx2)/(dx2*dy2)/180.0-(21*dy2+200*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(5,1) = -(140*dy2-36*dx2)/(dx2*dy2)/45.0-(20*dy2+24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(5,2) = -(140*dy2-36*dx2)/(dx2*dy2)/45.0-(20*dy2+24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(5,3) = -(40*dy2+36*dx2)/(dx2*dy2)/45.0-(40*dy2-24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(5,4) = -(40*dy2+36*dx2)/(dx2*dy2)/45.0-(40*dy2-24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(5,5) = 2.0*((140*dy2+24*dx2)/3.0+(20*dy2+24*dx2)/3.0)/(15.0*dx2*dy2);
    refDNjDNk(5,6) = 0;
    refDNjDNk(5,7) = 4.0*(40*dy2-24*dx2)/(45.0*dx2*dy2);
    refDNjDNk(5,8) = 0;
    refDNjDNk(6,1) = (21*dy2-160*dx2)/(dx2*dy2)/180.0-(69*dy2+160*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(6,2) = (69*dy2-440*dx2)/(dx2*dy2)/180.0-(21*dy2+200*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(6,3) = (69*dy2-440*dx2)/(dx2*dy2)/180.0-(21*dy2+200*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(6,4) = (21*dy2-160*dx2)/(dx2*dy2)/180.0-(69*dy2+160*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(6,5) = 0;
    refDNjDNk(6,6) = 2.0*(48*dy2+160*dx2)/(45.0*dx2*dy2);
    refDNjDNk(6,7) = 0;
    refDNjDNk(6,8) = (-4.0)*(24*dy2-40*dx2)/(45.0*dx2*dy2);
    refDNjDNk(7,1) = -(40*dy2+36*dx2)/(dx2*dy2)/45.0-(40*dy2-24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(7,2) = -(40*dy2+36*dx2)/(dx2*dy2)/45.0-(40*dy2-24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(7,3) = -(140*dy2-36*dx2)/(dx2*dy2)/45.0-(20*dy2+24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(7,4) = -(140*dy2-36*dx2)/(dx2*dy2)/45.0-(20*dy2+24*dx2)/(dx2*dy2)/45.0;
    refDNjDNk(7,5) = 4.0*(40*dy2-24*dx2)/(45.0*dx2*dy2);
    refDNjDNk(7,6) = 0;
    refDNjDNk(7,7) = 2.0*((140*dy2+24*dx2)/3.0+(20*dy2+24*dx2)/3.0)/(15.0*dx2*dy2);
    refDNjDNk(7,8) = 0;
    refDNjDNk(8,1) = (69*dy2-440*dx2)/(dx2*dy2)/180.0-(21*dy2+200*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(8,2) = (21*dy2-160*dx2)/(dx2*dy2)/180.0-(69*dy2+160*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(8,3) = (21*dy2-160*dx2)/(dx2*dy2)/180.0-(69*dy2+160*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(8,4) = (69*dy2-440*dx2)/(dx2*dy2)/180.0-(21*dy2+200*dx2)/(dx2*dy2)/180.0;
    refDNjDNk(8,5) = 0;
    refDNjDNk(8,6) = (-4.0)*(24*dy2-40*dx2)/(45.0*dx2*dy2);
    refDNjDNk(8,7) = 0;
    refDNjDNk(8,8) = 2.0*(48*dy2+160*dx2)/(45.0*dx2*dy2);

// scale to bring this into physical space
    refDNjDNk *= 0.5*dx*0.5*dy;
  }
}
