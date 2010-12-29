/**
 * @file	lcrectcoordsys.cxx
 *
 * @brief	Unit tests for Lucee::Rectcoordsys class
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// lucee includes
#include <LcAlignedRectCoordSys.h>
#include <LcTest.h>

void
test_1()
{
  Lucee::Vec3 xu, yu, zu;
  Lucee::AlignedRectCoordSys tx(0);
  tx.fillWithUnitVecs(xu, yu, zu);

  LC_ASSERT("Testing unit vector", xu[0] == 1);
  LC_ASSERT("Testing if xu is a unit vector", epsCmp(xu.getNorm(), 1.0));
  LC_ASSERT("Testing if yu is a unit vector", epsCmp(yu.getNorm(), 1.0));
  LC_ASSERT("Testing if zu is a unit vector", epsCmp(zu.getNorm(), 1.0));

  Lucee::Vec3 yuxzu = yu.cross(zu);
  for (unsigned i=0; i<3; ++i)
    LC_ASSERT("Testing righthanded-ness", yuxzu[i] == xu[i]);

  Lucee::AlignedRectCoordSys ty(1);
  ty.fillWithUnitVecs(xu, yu, zu);

  LC_ASSERT("Testing unit vector", xu[1] == 1);
  LC_ASSERT("Testing if xu is a unit vector", epsCmp(xu.getNorm(), 1.0));
  LC_ASSERT("Testing if yu is a unit vector", epsCmp(yu.getNorm(), 1.0));
  LC_ASSERT("Testing if zu is a unit vector", epsCmp(zu.getNorm(), 1.0));

  yuxzu = yu.cross(zu);
  for (unsigned i=0; i<3; ++i)
    LC_ASSERT("Testing righthanded-ness", yuxzu[i] == xu[i]);

  Lucee::AlignedRectCoordSys tz(2);
  tz.fillWithUnitVecs(xu, yu, zu);

  LC_ASSERT("Testing unit vector", xu[2] == 1);
  LC_ASSERT("Testing if xu is a unit vector", epsCmp(xu.getNorm(), 1.0));
  LC_ASSERT("Testing if yu is a unit vector", epsCmp(yu.getNorm(), 1.0));
  LC_ASSERT("Testing if zu is a unit vector", epsCmp(zu.getNorm(), 1.0));

  yuxzu = yu.cross(zu);
  for (unsigned i=0; i<3; ++i)
    LC_ASSERT("Testing righthanded-ness", yuxzu[i] == xu[i]);
}

void
test_2()
{
  Lucee::AlignedRectCoordSys tx(0);
  double u[3] = {1.0, 2.0, 3.0}, v[3], uu[3];
  
  tx.rotateVecToLocal(u, v);
  LC_ASSERT("Testing rotation of vector", v[0] == 1.0);
  LC_ASSERT("Testing rotation of vector", v[1] == 2.0);
  LC_ASSERT("Testing rotation of vector", v[2] == 3.0);

  tx.rotateVecToGlobal(v, uu);
  LC_ASSERT("Testing tx rotation back to global", arraycmp(u, uu, 3));

  Lucee::AlignedRectCoordSys ty(1);
  ty.rotateVecToLocal(u, v);
  LC_ASSERT("Testing rotation of vector", v[0] == 2.0);
  LC_ASSERT("Testing rotation of vector", v[1] == -1.0);
  LC_ASSERT("Testing rotation of vector", v[2] == 3.0);

  ty.rotateVecToGlobal(v, uu);
  LC_ASSERT("Testing ty rotation back to global", arraycmp(u, uu, 3));

  Lucee::AlignedRectCoordSys tz(2);
  tz.rotateVecToLocal(u, v);
  LC_ASSERT("Testing rotation of vector", v[0] == 3.0);
  LC_ASSERT("Testing rotation of vector", v[1] == 2.0);
  LC_ASSERT("Testing rotation of vector", v[2] == -1.0);

  tz.rotateVecToGlobal(v, uu);
  LC_ASSERT("Testing tz rotation back to global", arraycmp(u, uu, 3));
}

int
main(void) 
{
  LC_BEGIN_TESTS("lcrectcoordsys");
  test_1();
  test_2();
  LC_END_TESTS;
}
