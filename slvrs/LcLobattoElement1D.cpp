/**
 * @file	LcLobattoElement1D.cpp
 *
 * @brief       Reference finite element with serendipity basis
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLobattoElement1D.h>
#include <LcStructuredGridBase.h>

// std includes
#include <cmath>

namespace Lucee
{
  const char *LobattoElement1D::id = "Lobatto1D";

  LobattoElement1D::LobattoElement1D()
    : Lucee::NodalFiniteElementIfc(2), polyOrder(1), 
      refNjNk(2,2), refDNjDNk(2,2)
  {
  }

  void
  LobattoElement1D::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::NodalFiniteElementIfc::readInput(tbl);

// read in polynomial order
    polyOrder = tbl.getNumber("polyOrder");

    if (polyOrder == 1)
      this->setNumNodes(2);
    else if (polyOrder == 2)
      this->setNumNodes(3);
    else if (polyOrder == 3)
      this->setNumNodes(4);
    else
    {
      Lucee::Except lce("LobattoElement1D: Order must be 1, 2 or 3.");
      lce << " Provided " << polyOrder << " instead";
      throw lce;
    }

// initialize matrices
    if (polyOrder == 1)
      setupPoly1();
    else if (polyOrder == 2)
      setupPoly2();
    else if (polyOrder == 3)
      setupPoly3();
  }

  void
  LobattoElement1D::getLocalToGlobal(std::vector<int>& lgMap) const
  {
    int ix = this->currIdx[0];
    if (polyOrder == 1)
    { // two nodes per cell
      lgMap[0] = ix;
      lgMap[1] = ix+1;
    }
    else if (polyOrder == 2)
    { // three nodes per cell
      lgMap[0] = 2*ix;
      lgMap[1] = 2*ix+1;
      lgMap[2] = 2*ix+2;
    }
    else if (polyOrder == 3)
    { // four nodes per cell
      lgMap[0] = 3*ix;
      lgMap[1] = 3*ix+1;
      lgMap[2] = 3*ix+2;
      lgMap[3] = 3*ix+3;
    }
  }

  void
  LobattoElement1D::getMassMatrix(Lucee::Matrix<double> NjNk) const
  {
    NjNk = refNjNk;
  }

  void
  LobattoElement1D::getStiffnessMatrix(Lucee::Matrix<double> DNjDNk) const
  {
    DNjDNk = refDNjDNk;
  }

  void
  LobattoElement1D::setupPoly1()
  {
    unsigned shape[2] = {2,2};
    int start[2] = {1,1};

// get hold of grid
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();

// get grid spacing (this is assumed to be uniform for now)
    double dx = grid.getDx(0);

// mass matrix (automatically generated. See scripts/nodal-basis-1d.mac)
    refNjNk = Lucee::Matrix<double>(shape, start);
    refNjNk(1,1) = 2.0/3.0;
    refNjNk(1,2) = 1.0/3.0;
    refNjNk(2,1) = 1.0/3.0;
    refNjNk(2,2) = 2.0/3.0;

// scale to bring this into physical space
    refNjNk *= 0.5*dx;

// stiffness matrix (automatically generated. See scripts/nodal-basis-1d.mac)
    refDNjDNk = Lucee::Matrix<double>(shape, start);
    refDNjDNk(1,1) = 1.0/2.0;
    refDNjDNk(1,2) = (-1.0)/2.0;
    refDNjDNk(2,1) = (-1.0)/2.0;
    refDNjDNk(2,2) = 1.0/2.0;

// scale to bring this into physical space
    refDNjDNk *= 2/dx;
  }

  void
  LobattoElement1D::setupPoly2()
  {
    unsigned shape[2] = {3,3};
    int start[2] = {1,1};

// get hold of grid
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();

// get grid spacing (this is assumed to be uniform for now)
    double dx = grid.getDx(0);

// mass matrix (automatically generated. See scripts/nodal-basis-1d.mac)
    refNjNk = Lucee::Matrix<double>(shape, start);
    refNjNk(1,1) = 4.0/15.0;
    refNjNk(1,2) = 2.0/15.0;
    refNjNk(1,3) = (-1.0)/15.0;
    refNjNk(2,1) = 2.0/15.0;
    refNjNk(2,2) = 16.0/15.0;
    refNjNk(2,3) = 2.0/15.0;
    refNjNk(3,1) = (-1.0)/15.0;
    refNjNk(3,2) = 2.0/15.0;
    refNjNk(3,3) = 4.0/15.0;

// scale to bring this into physical space
    refNjNk *= 0.5*dx;

// stiffness matrix (automatically generated. See scripts/nodal-basis-1d.mac)
    refDNjDNk = Lucee::Matrix<double>(shape, start);
    refDNjDNk(1,1) = 7.0/6.0;
    refDNjDNk(1,2) = (-4.0)/3.0;
    refDNjDNk(1,3) = 1.0/6.0;
    refDNjDNk(2,1) = (-4.0)/3.0;
    refDNjDNk(2,2) = 8.0/3.0;
    refDNjDNk(2,3) = (-4.0)/3.0;
    refDNjDNk(3,1) = 1.0/6.0;
    refDNjDNk(3,2) = (-4.0)/3.0;
    refDNjDNk(3,3) = 7.0/6.0;

// scale to bring this into physical space
    refDNjDNk *= 2/dx;
  }

  void
  LobattoElement1D::setupPoly3()
  {
    unsigned shape[2] = {4,4};
    int start[2] = {1,1};

// get hold of grid
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();

// get grid spacing (this is assumed to be uniform for now)
    double dx = grid.getDx(0);

// mass matrix (automatically generated. See scripts/nodal-basis-1d.mac)
    refNjNk = Lucee::Matrix<double>(shape, start);
    refNjNk(1,1) = 1.0/7.0;
    refNjNk(1,2) = std::sqrt(5)/42.0;
    refNjNk(1,3) = -std::sqrt(5)/42.0;
    refNjNk(1,4) = 1.0/42.0;
    refNjNk(2,1) = std::sqrt(5)/42.0;
    refNjNk(2,2) = 5.0/7.0;
    refNjNk(2,3) = 5.0/42.0;
    refNjNk(2,4) = -std::sqrt(5)/42.0;
    refNjNk(3,1) = -std::sqrt(5)/42.0;
    refNjNk(3,2) = 5.0/42.0;
    refNjNk(3,3) = 5.0/7.0;
    refNjNk(3,4) = std::sqrt(5)/42.0;
    refNjNk(4,1) = 1.0/42.0;
    refNjNk(4,2) = -std::sqrt(5)/42.0;
    refNjNk(4,3) = std::sqrt(5)/42.0;
    refNjNk(4,4) = 1.0/7.0;

// scale to bring this into physical space
    refNjNk *= 0.5*dx;

// stiffness matrix (automatically generated. See scripts/nodal-basis-1d.mac)
    refDNjDNk = Lucee::Matrix<double>(shape, start);
    refDNjDNk(1,1) = 13.0/6.0;
    refDNjDNk(1,2) = -(3*std::pow(5.0, 1.5)+25)/24.0;
    refDNjDNk(1,3) = (3*std::pow(5.0, 1.5)-25)/24.0;
    refDNjDNk(1,4) = (-1.0)/12.0;
    refDNjDNk(2,1) = -(3*std::pow(5.0, 1.5)+25)/24.0;
    refDNjDNk(2,2) = 25.0/6.0;
    refDNjDNk(2,3) = (-25.0)/12.0;
    refDNjDNk(2,4) = (3*std::pow(5.0, 1.5)-25)/24.0;
    refDNjDNk(3,1) = (3*std::pow(5.0, 1.5)-25)/24.0;
    refDNjDNk(3,2) = (-25.0)/12.0;
    refDNjDNk(3,3) = 25.0/6.0;
    refDNjDNk(3,4) = -(3*std::pow(5.0, 1.5)+25)/24.0;
    refDNjDNk(4,1) = (-1.0)/12.0;
    refDNjDNk(4,2) = (3*std::pow(5.0, 1.5)-25)/24.0;
    refDNjDNk(4,3) = -(3*std::pow(5.0, 1.5)+25)/24.0;
    refDNjDNk(4,4) = 13.0/6.0;

// scale to bring this into physical space
    refDNjDNk *= 2/dx;
  }
}
