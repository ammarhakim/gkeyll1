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
#include <LcGaussianQuadRule.h>
#include <LcLobattoElement1D.h>
#include <LcStructuredGridBase.h>

// std includes
#include <cmath>

namespace Lucee
{
  const char *LobattoElement1D::id = "Lobatto";

  static
  double
  mapUnitToX(double xc, double dx, double eta)
  {
    return 0.5*dx*eta + xc;
  }

  LobattoElement1D::LobattoElement1D()
    : Lucee::NodalFiniteElementIfc<1>(2), polyOrder(1), 
      refNjNk(2,2), refDNjDNk(2,2), refDNjNk_0(2,2)
  {
  }

  void
  LobattoElement1D::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::NodalFiniteElementIfc<1>::readInput(tbl);

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

// determine number of global degrees of freedom
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();
    unsigned gvol = grid.getGlobalRegion().getVolume();

    if (polyOrder == 1)
      numGlobal = gvol+1; // 1 *global* node per cell and one extra for right-most node
    else if (polyOrder == 2)
      numGlobal = 2*gvol+1; // 2 *global* nodes per cell and one extra for right-mode node
    else if (polyOrder == 3)
      numGlobal = 3*gvol+1; // 3 *global* nodes per cell and one extra for right-mode node
  }

  unsigned
  LobattoElement1D::getNumSurfLowerNodes(unsigned dir) const
  {
    return 1;
  }

  unsigned
  LobattoElement1D::getNumSurfUpperNodes(unsigned dir) const
  {
    return 1;
  }

  unsigned
  LobattoElement1D::getNumGlobalNodes() const
  {
    return numGlobal;
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
  LobattoElement1D::getSurfLowerLocalToGlobal(unsigned dir,
    std::vector<int>& lgMap) const
  {
    int ix = this->currIdx[0];
    if (polyOrder == 1)
    { // two nodes per cell
      lgMap[0] = ix;
    }
    else if (polyOrder == 2)
    { // three nodes per cell
      lgMap[0] = 2*ix;
    }
    else if (polyOrder == 3)
    { // four nodes per cell
      lgMap[0] = 3*ix;
    }
  }

  void
  LobattoElement1D::getSurfUpperLocalToGlobal(unsigned dir,
    std::vector<int>& lgMap) const
  {
    int ix = this->currIdx[0];
    if (polyOrder == 1)
    { // two nodes per cell
      lgMap[0] = ix+1;
    }
    else if (polyOrder == 2)
    { // three nodes per cell
      lgMap[0] = 2*ix+2;
    }
    else if (polyOrder == 3)
    { // four nodes per cell
      lgMap[0] = 3*ix+3;
    }
  }

  void
  LobattoElement1D::getSurfLowerNodeNums(unsigned dir, std::vector<int>& nodeNum) const
  {
    nodeNum[0] = 0;
  }

  void
  LobattoElement1D::getSurfUpperNodeNums(unsigned dir, std::vector<int>& nodeNum) const
  {
    nodeNum[0] = this->getNumNodes()-1;
  }

  void
  LobattoElement1D::getExclusiveNodeIndices(std::vector<int>& ndIds)
  {
    ndIds.clear();
    ndIds.resize(this->getNumNodes()-1);
    for (unsigned n=0; n<this->getNumNodes()-1; ++n)
      ndIds[n] = n;
  }

  void
  LobattoElement1D::getNodalCoordinates(Lucee::Matrix<double>& nodeCoords)
  {
// get grid
    const Lucee::StructuredGridBase<1>& grid
      = this->getGrid<Lucee::StructuredGridBase<1> >();
// set index and get centroid coordinate
    grid.setIndex(this->currIdx);
    double xc[3], dx;
    grid.getCentroid(xc);
    dx = grid.getDx(0);

    nodeCoords = 0.0;
    if (polyOrder == 1)
    {
      nodeCoords(0,0) = mapUnitToX(xc[0], dx, -1);
      nodeCoords(1,0) = mapUnitToX(xc[0], dx, 1);
    }
    else if (polyOrder == 2)
    {
      nodeCoords(0,0) = mapUnitToX(xc[0], dx, -1);
      nodeCoords(1,0) = mapUnitToX(xc[0], dx, 0);
      nodeCoords(2,0) = mapUnitToX(xc[0], dx, 1);
    }
    else if (polyOrder == 3)
    {
      double s5 = std::sqrt(1.0/5.0);
      nodeCoords(0,0) = mapUnitToX(xc[0], dx, -1);
      nodeCoords(1,0) = mapUnitToX(xc[0], dx, -s5);
      nodeCoords(2,0) = mapUnitToX(xc[0], dx, s5);
      nodeCoords(3,0) = mapUnitToX(xc[0], dx, 1);
    }
  }

  void
  LobattoElement1D::getWeights(std::vector<double>& w)
  {
    unsigned nn = this->getNumNodes();
    for (unsigned k=0; k<nn; ++k) w[k] = weights[k];
  }

  void
  LobattoElement1D::getMassMatrix(Lucee::Matrix<double>& NjNk) const
  {
    NjNk.copy(refNjNk);
  }

  void
  LobattoElement1D::getStiffnessMatrix(Lucee::Matrix<double>& DNjDNk) const
  {
    DNjDNk.copy(refDNjDNk);
  }

  void
  LobattoElement1D::getGradStiffnessMatrix(unsigned dir, Lucee::Matrix<double>& DNjNk) const
  {
    if (dir>0)
      throw Lucee::Except(
        "getGradStiffnessMatrix::getGradStiffnessMatrix: Can only use this basis in 1D!");
    DNjNk.copy(refDNjNk_0);
  }

  unsigned
  LobattoElement1D::getNumGaussNodes() const
  {
    return polyOrder+1;
  }

  void
  LobattoElement1D::getGaussQuadData(Lucee::Matrix<double>& interpMat,
    Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const
  {
// copy over data
    interpMat.copy(gaussData.interpMat);
    ordinates.copy(gaussData.ords);
    for (unsigned i=0; i<gaussData.weights.size(); ++i)
      weights[i] = gaussData.weights[i];
  }

  void
  LobattoElement1D::extractFromField(const Lucee::Field<1, double>& fld,
    std::vector<double>& data)
  {
    Lucee::ConstFieldPtr<double> fldPtr = fld.createConstPtr();
    Lucee::ConstFieldPtr<double> fldPtrp = fld.createConstPtr();
// attach pointers to proper locations
    fld.setPtr(fldPtr, this->currIdx[0]);
    fld.setPtr(fldPtrp, this->currIdx[0]+1);

    unsigned nlocal = this->getNumNodes();
// extract data
      for (unsigned k=0; k<nlocal-1; ++k)
        data[k] = fldPtr[k];
      data[nlocal-1] = fldPtrp[0]; // get right-most data from first node of right cell
  }

  void
  LobattoElement1D::copyAllDataFromField(const Lucee::Field<1, double>& fld,
    double *data)
  {
// region to copy
    Lucee::Region<1, int> rgn =
      this->getGrid<Lucee::StructuredGridBase<1> >().getLocalRegion();

    unsigned nlocal = this->getNumNodes();
    Lucee::ConstFieldPtr<double> fldPtr = fld.createConstPtr();
// copy data
    unsigned count = 0;
    for (int i=rgn.getLower(0); i<rgn.getUpper(0); ++i)
    {
      fld.setPtr(fldPtr, i);
      for (unsigned k=0; k<nlocal-1; ++k)
        data[count++] = fldPtr[k];
    }
// copy data at last node
    fld.setPtr(fldPtr, rgn.getUpper(0));
    data[count] = fldPtr[0];
  }

  void
  LobattoElement1D::copyAllDataToField(const double *data, Lucee::Field<1, double>& fld)
  {
// region to copy
    Lucee::Region<1, int> rgn =
      this->getGrid<Lucee::StructuredGridBase<1> >().getLocalRegion();

    unsigned nlocal = this->getNumNodes();
    Lucee::FieldPtr<double> fldPtr = fld.createPtr();
// copy data
    unsigned count = 0;
    for (int i=rgn.getLower(0); i<rgn.getUpper(0); ++i)
    {
      fld.setPtr(fldPtr, i);
      for (unsigned k=0; k<nlocal-1; ++k)
        fldPtr[k] = data[count++];
    }
// copy field at last node
    fld.setPtr(fldPtr, rgn.getUpper(0));
    fldPtr[0] = data[count];
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

// grad-stiffness matrix (automatically generated. See scripts/nodal-basis-1d.mac)
    refDNjNk_0 = Lucee::Matrix<double>(shape, start);
    refDNjNk_0(1,1) = -1/dx;
    refDNjNk_0(1,2) = -1/dx;
    refDNjNk_0(2,1) = 1/dx;
    refDNjNk_0(2,2) = 1/dx;

// scale to bring this into physical space
    refDNjNk_0 *= 0.5*dx;

// compute weights
    weights.resize(2);
    weights[0] = 1.0;
    weights[1] = 1.0;

    for (unsigned i=0; i<2; ++i)
      weights[i] = 0.5*dx*weights[i];

// Gaussian quadrature data
    calcGuassData(2);

// interpolation matrix (automatically generated. See scripts/nodal-basis-1d.mac)
    gaussData.interpMat(1,1) = 1/std::sqrt(3.0)/2.0+1.0/2.0;
    gaussData.interpMat(1,2) = 1.0/2.0-1/std::sqrt(3.0)/2.0;
    gaussData.interpMat(2,1) = 1.0/2.0-1/std::sqrt(3.0)/2.0;
    gaussData.interpMat(2,2) = 1/std::sqrt(3.0)/2.0+1.0/2.0;
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

// grad-stiffness matrix (automatically generated. See scripts/nodal-basis-1d.mac)
    refDNjNk_0 = Lucee::Matrix<double>(shape, start);
    refDNjNk_0(1,1) = -1/dx;
    refDNjNk_0(1,2) = (-4.0)/(3.0*dx);
    refDNjNk_0(1,3) = 1/dx/3.0;
    refDNjNk_0(2,1) = 4.0/(3.0*dx);
    refDNjNk_0(2,2) = 0;
    refDNjNk_0(2,3) = (-4.0)/(3.0*dx);
    refDNjNk_0(3,1) = -1/dx/3.0;
    refDNjNk_0(3,2) = 4.0/(3.0*dx);
    refDNjNk_0(3,3) = 1/dx;

// scale to bring this into physical space
    refDNjNk_0 *= 0.5*dx;

// compute weights
    weights.resize(3);
    weights[0] = 1.0/3.0;
    weights[1] = 4.0/3.0;
    weights[2] = 1.0/3.0;

    for (unsigned i=0; i<3; ++i)
      weights[i] = 0.5*dx*weights[i];

// Gaussian quadrature data
    calcGuassData(3);

// interpolation matrix (automatically generated. See scripts/nodal-basis-1d.mac)
    gaussData.interpMat(1,1) = .6872983346207417;
    gaussData.interpMat(1,2) = .3999999999999997;
    gaussData.interpMat(1,3) = -.08729833462074155;
    gaussData.interpMat(2,1) = 0.0;
    gaussData.interpMat(2,2) = 1.0;
    gaussData.interpMat(2,3) = 0.0;
    gaussData.interpMat(3,1) = -.08729833462074174;
    gaussData.interpMat(3,2) = .4000000000000001;
    gaussData.interpMat(3,3) = .6872983346207415;
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
    refNjNk(1,2) = std::sqrt(5.0)/42.0;
    refNjNk(1,3) = -std::sqrt(5.0)/42.0;
    refNjNk(1,4) = 1.0/42.0;
    refNjNk(2,1) = std::sqrt(5.0)/42.0;
    refNjNk(2,2) = 5.0/7.0;
    refNjNk(2,3) = 5.0/42.0;
    refNjNk(2,4) = -std::sqrt(5.0)/42.0;
    refNjNk(3,1) = -std::sqrt(5.0)/42.0;
    refNjNk(3,2) = 5.0/42.0;
    refNjNk(3,3) = 5.0/7.0;
    refNjNk(3,4) = std::sqrt(5.0)/42.0;
    refNjNk(4,1) = 1.0/42.0;
    refNjNk(4,2) = -std::sqrt(5.0)/42.0;
    refNjNk(4,3) = std::sqrt(5.0)/42.0;
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

// grad-stiffness matrix (automatically generated. See scripts/nodal-basis-1d.mac)
    refDNjNk_0 = Lucee::Matrix<double>(shape, start);

    refDNjNk_0(1,1) = -1/dx;
    refDNjNk_0(1,2) = -(std::pow(5, 1.5)+5)/dx/12.0;
    refDNjNk_0(1,3) = (std::pow(5, 1.5)-5)/dx/12.0;
    refDNjNk_0(1,4) = -1/dx/6.0;
    refDNjNk_0(2,1) = (std::pow(5, 1.5)+5)/dx/12.0;
    refDNjNk_0(2,2) = 0;
    refDNjNk_0(2,3) = -std::pow(5, 1.5)/dx/6.0;
    refDNjNk_0(2,4) = (std::pow(5, 1.5)-5)/dx/12.0;
    refDNjNk_0(3,1) = -(std::pow(5, 1.5)-5)/dx/12.0;
    refDNjNk_0(3,2) = std::pow(5, 1.5)/dx/6.0;
    refDNjNk_0(3,3) = 0;
    refDNjNk_0(3,4) = -(std::pow(5, 1.5)+5)/dx/12.0;
    refDNjNk_0(4,1) = 1/dx/6.0;
    refDNjNk_0(4,2) = -(std::pow(5, 1.5)-5)/dx/12.0;
    refDNjNk_0(4,3) = (std::pow(5, 1.5)+5)/dx/12.0;
    refDNjNk_0(4,4) = 1/dx;

// scale to bring this into physical space
    refDNjNk_0 *= 0.5*dx;

// compute weights
    weights.resize(4);
    double c5 = std::pow(5.0, 1.5);

    weights[0] = 1.0/6.0;
    weights[1] = (3*c5+40)/96.0-(3*c5-40)/96.0;
    weights[2] = (3*c5+40)/96.0-(3*c5-40)/96.0;
    weights[3] = 1.0/6.0;

    for (unsigned i=0; i<2; ++i)
      weights[i] = 0.5*dx*weights[i];

// Gaussian quadrature data
    calcGuassData(4);

// interpolation matrix (automatically generated. See scripts/nodal-basis-1d.mac)
    gaussData.interpMat(1,1) = .6299431661034449;
    gaussData.interpMat(1,2) = .4725587471138189;
    gaussData.interpMat(1,3) = -0.14950343104608;
    gaussData.interpMat(1,4) = .04700151782881631;
    gaussData.interpMat(2,1) = -.07069479527385608;
    gaussData.interpMat(2,2) = .9729761862582639;
    gaussData.interpMat(2,3) = 0.132539926245426;
    gaussData.interpMat(2,4) = -.03482131722983343;
    gaussData.interpMat(3,1) = -.03482131722983417;
    gaussData.interpMat(3,2) = .1325399262454268;
    gaussData.interpMat(3,3) = .9729761862582633;
    gaussData.interpMat(3,4) = -.07069479527385587;
    gaussData.interpMat(4,1) = 0.0470015178288162;
    gaussData.interpMat(4,2) = -.1495034310460798;
    gaussData.interpMat(4,3) = .4725587471138186;
    gaussData.interpMat(4,4) = .6299431661034451;
  }

  void
  LobattoElement1D::calcGuassData(unsigned nord)
  {
// get hold of grid
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();

// get grid spacing (this is assumed to be uniform for now)
    double dx[1];
    dx[0] = grid.getDx(0);

    Lucee::GaussianQuadRule gauss(nord);
    Lucee::Vector<double> w(nord), mu(nord);
    gauss.getOrdinatesAndWeights(mu, w);

// allocate memory for quadrature
    gaussData.reset(nord, nord);

// set ordinates
    for (unsigned i=0; i<nord; ++i)
    {
      gaussData.ords(i,0) = mu[i];
      gaussData.ords(i,1) = 0.0;
      gaussData.ords(i,2) = 0.0;
    }
   
// set weights
    for (unsigned i=0; i<nord; ++i)
      gaussData.weights[i] = 0.5*dx[0]*w[i];
  }
}
