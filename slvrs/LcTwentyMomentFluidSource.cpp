/**
 * @file	LcTwentyMomentFluidSource.cpp
 *
 * @brief       Compute source terms in pressure tensor equations.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcTwentyMomentFluidSource.h>

namespace Lucee
{
  // makes indexing easier
  static const unsigned RHO = 0;
  static const unsigned U1 = 1;
  static const unsigned U2 = 2;
  static const unsigned U3 = 3;
  static const unsigned P11 = 4;
  static const unsigned P12 = 5;
  static const unsigned P13 = 6;
  static const unsigned P22 = 7;
  static const unsigned P23 = 8;
  static const unsigned P33 = 9;
  static const unsigned Q111 = 10;
  static const unsigned Q112 = 11;
  static const unsigned Q113 = 12;
  static const unsigned Q122 = 13;
  static const unsigned Q123 = 14;
  static const unsigned Q133 = 15;
  static const unsigned Q222 = 16;
  static const unsigned Q223 = 17;
  static const unsigned Q233 = 18;
  static const unsigned Q333 = 19;

  // set id
  const char *TwentyMomentFluidSource::id = "TwentyMomentFluid";
  
  TwentyMomentFluidSource::TwentyMomentFluidSource()
    : Lucee::PointSourceIfc(26,20)
  {
  }

  void TwentyMomentFluidSource::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::PointSourceIfc::readInput(tbl);

    double charge = tbl.getNumber("charge");
    double mass = tbl.getNumber("mass");
    qbym = charge/mass;

    hasCollisions = false;
    if (tbl.hasNumber("collisionFrequency"))
    {
      hasCollisions = true;
      nu = tbl.getNumber("collisionFrequency");
    }
  }

  void TwentyMomentFluidSource::getSource(double tm, const double loc [3], std::vector<double>& src)
  {
    double r = this->getData(0);
    double u = this->getData(1)/r;
    double v = this->getData(2)/r;
    double w = this->getData(3)/r;
    double pxx = this->getData(4);
    double pxy = this->getData(5);
    double pxz = this->getData(6);
    double pyy = this->getData(7);
    double pyz = this->getData(8);
    double pzz = this->getData(9);
    double qxxx = this->getData(10);
    double qxxy = this->getData(11);
    double qxxz = this->getData(12);
    double qxyy = this->getData(13);
    double qxyz = this->getData(14);
    double qxzz = this->getData(15);
    double qyyy = this->getData(16);
    double qyyz = this->getData(17);
    double qyzz = this->getData(18);
    double qzzz = this->getData(19);
    double ex = this->getData(20);
    double ey = this->getData(21);
    double ez = this->getData(22);
    double bx = this->getData(23);
    double by = this->getData(24);
    double bz = this->getData(25);

    double re = qbym;

    src[RHO] = 0.0;
// momentum source terms
    src[U1] = r*re*(ex+v*bz-w*by);
    src[U2] = r*re*(ey+w*bx-u*bz);
    src[U3] = r*re*(ez+u*by-v*bx);
// pressure tensor source terms
    src[P11] = 2.0*r*re*u*ex+2.0*re*(bz*pxy-by*pxz);
    src[P12] = r*re*(u*ey+v*ex)+re*(bz*pyy-by*pyz-bz*pxx+bx*pxz);
    src[P13] = r*re*(u*ez+w*ex)+re*(bz*pyz+by*pxx-by*pzz-bx*pxy);
    src[P22] = 2.0*r*re*v*ey+2.0*re*(bx*pyz-bz*pxy);
    src[P23] = r*re*(v*ez+w*ey)+re*(by*pxy-bz*pxz+bx*pzz-bx*pyy);
    src[P33] = 2.0*r*re*w*ez+2.0*re*(by*pxz-bx*pyz);
// heat flux tensor source terms. 
// dQ/dt + div R = q/m sym(E_i Pjk) + - q/m sym(B cross Q)
    src[Q111] = re*(pxx*ex + pxx*ex + pxx*ex) - re*(by*qxxz-bz*qxxy + by*qxxz-bz*qxxy + by*qxxz-bz*qxxy);
    src[Q112] = re*(pxx*ey + pxy*ex + pxy*ex) - re*(by*qxyz-bz*qxyy + by*qxyz-bz*qxyy + bz*qxxx-bx*qxxz);
    src[Q113] = re*(pxx*ez + pxz*ex + pxz*ex) - re*(by*qxzz-bz*qxyz + by*qxzz-bz*qxyz + bx*qxxy-by*qxxx);
    src[Q122] = re*(pxy*ey + pxy*ey + pyy*ex) - re*(by*qyyz-bz*qyyy + bz*qxxy-bx*qxyz + bz*qxxy-bx*qxyz);
    src[Q123] = re*(pxy*ez + pxz*ey + pyz*ex) - re*(by*qyzz-bz*qyyz + bz*qxxz-bx*qxzz + bx*qxyy-by*qxxy);
    src[Q133] = re*(pxz*ez + pxz*ez + pzz*ex) - re*(by*qzzz-bz*qyzz + bx*qxyz-by*qxxz + bx*qxyz-by*qxxz);
    src[Q222] = re*(pyy*ey + pyy*ey + pyy*ey) - re*(bz*qxyy-bx*qyyz + bz*qxyy-bx*qyyz + bz*qxyy-bx*qyyz);
    src[Q223] = re*(pyy*ez + pyz*ey + pyz*ey) - re*(bz*qxyz-bx*qyzz + bz*qxyz-bx*qyzz + bx*qyyy-by*qxyy);
    src[Q233] = re*(pyz*ez + pyz*ez + pzz*ey) - re*(bz*qxzz-bx*qzzz + bx*qyyz-by*qxyz + bx*qyyz-by*qxyz);
    src[Q333] = re*(pzz*ez + pzz*ez + pzz*ez) - re*(bx*qyzz-by*qxzz + bx*qyzz-by*qxzz + bx*qyzz-by*qxzz);
// ignore collisions for now.

    if (hasCollisions)
    {
// compute pressure tensor in fluid frame
      double fpxx = pxx-r*u*u;
      double fpxy = pxy-r*u*v;
      double fpxz = pxz-r*u*w;
      double fpyy = pyy-r*v*v;
      double fpyz = pyz-r*v*w;
      double fpzz = pzz-r*w*w;
      
// compute scalar pressure
      double pr = (fpxx+fpyy+fpzz)/3.0;

// add in collision terms
      src[P11] += nu*(pr-fpxx);
      src[P12] += nu*(0-fpxy);
      src[P13] += nu*(0-fpxz);
      src[P22] += nu*(pr-fpyy);
      src[P23] += nu*(0-fpyz);
      src[P33] += nu*(pr-fpzz); 
      double q111 = qxxx - 3*fpxx*u - r*u*u*u;
      double q112 = qxxy - 2*fpxy*u - fpxx*v - r*u*u*v;
      double q113 = qxxz - 2*fpxz*u - fpxx*w - r*u*u*w;
      double q122 = qxyy - fpyy*u - 2*fpxy*v - r*u*v*v;
      double q123 = qxyz - fpyz*u - fpxz*v - fpxy*w - r*u*v*w;
      double q133 = qxzz - fpzz*u - 2*fpxz*w - r*u*w*w;
      double q222 = qyyy - 3*fpyy*v - r*v*v*v;
      double q223 = qyyz - 2*fpyz*v - fpyy*w - r*v*v*w;
      double q233 = qyzz - fpzz*v - 2*fpyz*w - r*v*w*w;
      double q333 = qzzz - 3*fpzz*w - r*w*w*w;
      
      src[Q111] += -nu*q111;
      src[Q112] += -nu*q112;
      src[Q113] += -nu*q113;
      src[Q122] += -nu*q122;
      src[Q123] += -nu*q123;
      src[Q133] += -nu*q133;
      src[Q222] += -nu*q222;
      src[Q223] += -nu*q223;
      src[Q233] += -nu*q233;
      src[Q333] += -nu*q333;
    }

  }

  void TwentyMomentFluidSource::getSourceJac(double tm, const double loc [3], Lucee::Matrix<double>& jac)
  {
    throw Lucee::Except("TenMomentFluidSource::getSourceJac: Not implemented!");
  }


}

