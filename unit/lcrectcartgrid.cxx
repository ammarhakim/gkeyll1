/**
 * @file	lcrectcartgrid.cxx
 *
 * @brief	Unit tests for Lucee::RectCartGrid class
 */

// lucee includes
#include <LcTest.h>
#include <LcRectCartGrid.h>
#include <LcRowMajorSequencer.h>
#include <LcVec3.h>

// std includes
#include <cmath>

void
test_1()
{
  int cells[3] = {1, 2, 4};
  Lucee::Region<3, int> localBox(cells);
  double lower[3] = {0.5, 0.5, 0.5};
  double upper[3] = {1.5, 2.5, 4.5};
  Lucee::Region<3, double> physBox(lower, upper);
  Lucee::RectCartGrid<3> grid3(localBox, localBox, physBox);

// cell spacing
  LC_ASSERT("Testing spacing in x-direction", grid3.getDx(0) == 1.0);
  LC_ASSERT("Testing spacing in y-direction", grid3.getDx(1) == 1.0);
  LC_ASSERT("Testing spacing in z-direction", grid3.getDx(2) == 1.0);

// check shape of physical space
  Lucee::Region<3, double> pbox = grid3.getComputationalSpace();
  for (unsigned i=0; i<3; ++i)
  {
    LC_ASSERT("Testing computational box lower coordinate", 
      pbox.getLower(i) == physBox.getLower(i));
    LC_ASSERT("Testing computational box upper coordinate", 
      pbox.getUpper(i) == physBox.getUpper(i));
    LC_ASSERT("Testing computational shape", 
      pbox.getShape(i) == physBox.getShape(i));
  }

// cell volume
  int idx[3];
  Lucee::RowMajorSequencer<3> seq(localBox);
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid3.setIndex(idx);
    LC_ASSERT("Testing cell volume", grid3.getVolume() == 1.0);
  }

  LC_ASSERT("Testing number of cells", grid3.getNumCells(0) == 1);
  LC_ASSERT("Testing number of cells", grid3.getNumCells(1) == 2);
  LC_ASSERT("Testing number of cells", grid3.getNumCells(2) == 4);

  seq.reset();
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid3.setIndex(idx);
    for (unsigned i=0; i<3; ++i)
      LC_ASSERT("Testing surface area volume", grid3.getSurfArea(i) == 1.0);
  }

  double xc[3];
  seq.reset();
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid3.setIndex(idx);

    grid3.getCentroid(xc);
    for (unsigned i=0; i<3; ++i)
      LC_ASSERT("Testing centroid coordinates", xc[i] == 0.5+1.0*(idx[i]+0.5));

    grid3.getVertex(xc);
    for (unsigned i=0; i<3; ++i)
      LC_ASSERT("Testing vertex coordinates", xc[i] == 0.5+1.0*idx[i]);
  }

  double norm[3], tan1[3], tan2[3];
  Lucee::Vec3<double> ux(1, 0, 0), uy(0, 1, 0), uz(0, 0, 1);

  seq.reset();
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid3.setIndex(idx);

    grid3.getSurfCoordSys(0, norm, tan1, tan2);
    Lucee::Vec3<double> lx(tan1), ly(tan2), lz(norm), t1t2;
    
    LC_ASSERT("Testing normal", ux.dot(lz) == 1.0);
    t1t2 = lx.cross(ly);
    for (unsigned i=0; i<3; ++i)
      LC_ASSERT("Testing right-handedness of coord-sys", t1t2[i] == lz[i]);
  }

  seq.reset();
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid3.setIndex(idx);

    grid3.getSurfCoordSys(1, norm, tan1, tan2);
    Lucee::Vec3<double> lx(tan1), ly(tan2), lz(norm), t1t2;
    
    LC_ASSERT("Testing normal", uy.dot(lz) == 1.0);
    t1t2 = lx.cross(ly);
    for (unsigned i=0; i<3; ++i)
      LC_ASSERT("Testing right-handedness of coord-sys", t1t2[i] == lz[i]);
  }

  seq.reset();
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid3.setIndex(idx);

    grid3.getSurfCoordSys(2, norm, tan1, tan2);
    Lucee::Vec3<double> lx(tan1), ly(tan2), lz(norm), t1t2;
    
    LC_ASSERT("Testing normal", uz.dot(lz) == 1.0);
    t1t2 = lx.cross(ly);
    for (unsigned i=0; i<3; ++i)
      LC_ASSERT("Testing right-handedness of coord-sys", t1t2[i] == lz[i]);
  }

// write grid to file
  grid3.write("lcrectcartgrid-grid3.h5");

}

void
test_2()
{
  int cells[3] = {1, 4, 16};
  Lucee::Region<3, int> localBox(cells);
  double lower[3] = {0.0, 0.0, 0.0};
  double upper[3] = {1.0, 2.0, 4.0};
  Lucee::Region<3, double> physBox(lower, upper);
  Lucee::RectCartGrid<3> grid3(localBox, localBox, physBox);

// cell volume
  int idx[3];
  Lucee::RowMajorSequencer<3> seq(localBox);
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid3.setIndex(idx);
    LC_ASSERT("Testing cell volume", grid3.getVolume() == 1*0.5*0.25);
  }

  LC_ASSERT("Testing number of cells", grid3.getNumCells(0) == 1);
  LC_ASSERT("Testing number of cells", grid3.getNumCells(1) == 4);
  LC_ASSERT("Testing number of cells", grid3.getNumCells(2) == 16);

  seq.reset();
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid3.setIndex(idx);
    LC_ASSERT("Testing surface area volume", grid3.getSurfArea(0) == 0.5*0.25);
    LC_ASSERT("Testing surface area volume", grid3.getSurfArea(1) == 1.0*0.25);
    LC_ASSERT("Testing surface area volume", grid3.getSurfArea(2) == 1.0*0.5);
  }

  double xc[3];
  seq.reset();
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid3.setIndex(idx);

    grid3.getCentroid(xc);
    LC_ASSERT("Testing centroid coordinates", xc[0] == 1.0*(idx[0]+0.5));
    LC_ASSERT("Testing centroid coordinates", xc[1] == 0.5*(idx[1]+0.5));
    LC_ASSERT("Testing centroid coordinates", xc[2] == 0.25*(idx[2]+0.5));
  }

  double norm[3], tan1[3], tan2[3];
  Lucee::Vec3<double> ux(1, 0, 0), uy(0, 1, 0), uz(0, 0, 1);

  seq.reset();
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid3.setIndex(idx);

    grid3.getSurfCoordSys(0, norm, tan1, tan2);
    Lucee::Vec3<double> lx(tan1), ly(tan2), lz(norm), t1t2;
    
    LC_ASSERT("Testing normal", ux.dot(lz) == 1.0);
    t1t2 = lx.cross(ly);
    for (unsigned i=0; i<3; ++i)
      LC_ASSERT("Testing right-handedness of coord-sys", t1t2[i] == lz[i]);
  }

  seq.reset();
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid3.setIndex(idx);

    grid3.getSurfCoordSys(1, norm, tan1, tan2);
    Lucee::Vec3<double> lx(tan1), ly(tan2), lz(norm), t1t2;
    
    LC_ASSERT("Testing normal", uy.dot(lz) == 1.0);
    t1t2 = lx.cross(ly);
    for (unsigned i=0; i<3; ++i)
      LC_ASSERT("Testing right-handedness of coord-sys", t1t2[i] == lz[i]);
  }

  seq.reset();
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid3.setIndex(idx);

    grid3.getSurfCoordSys(2, norm, tan1, tan2);
    Lucee::Vec3<double> lx(tan1), ly(tan2), lz(norm), t1t2;
    
    LC_ASSERT("Testing normal", uz.dot(lz) == 1.0);
    t1t2 = lx.cross(ly);
    for (unsigned i=0; i<3; ++i)
      LC_ASSERT("Testing right-handedness of coord-sys", t1t2[i] == lz[i]);
  }

}

void
test_3()
{
  int cells[3] = {1, 4, 16};
  Lucee::Region<3, int> localBox(cells);
  double lower[3] = {-1.0, -2.0, -4.0};
  double upper[3] = {0.0, 0.0, 0.0};
  Lucee::Region<3, double> physBox(lower, upper);
  Lucee::RectCartGrid<3> grid3(localBox, localBox, physBox);

// cell volume
  int idx[3];
  Lucee::RowMajorSequencer<3> seq(localBox);
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid3.setIndex(idx);
    LC_ASSERT("Testing cell volume", grid3.getVolume() == 1*0.5*0.25);
  }

  LC_ASSERT("Testing number of cells", grid3.getNumCells(0) == 1);
  LC_ASSERT("Testing number of cells", grid3.getNumCells(1) == 4);
  LC_ASSERT("Testing number of cells", grid3.getNumCells(2) == 16);

  seq.reset();
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid3.setIndex(idx);
    LC_ASSERT("Testing surface area volume", grid3.getSurfArea(0) == 0.5*0.25);
    LC_ASSERT("Testing surface area volume", grid3.getSurfArea(1) == 1.0*0.25);
    LC_ASSERT("Testing surface area volume", grid3.getSurfArea(2) == 1.0*0.5);
  }

  double xc[3];
  seq.reset();
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid3.setIndex(idx);

    grid3.getCentroid(xc);
    LC_ASSERT("Testing centroid coordinates", xc[0] == -1.0 + 1.0*(idx[0]+0.5));
    LC_ASSERT("Testing centroid coordinates", xc[1] == -2.0 + 0.5*(idx[1]+0.5));
    LC_ASSERT("Testing centroid coordinates", xc[2] == -4.0 + 0.25*(idx[2]+0.5));
  }

  double norm[3], tan1[3], tan2[3];
  Lucee::Vec3<double> ux(1, 0, 0), uy(0, 1, 0), uz(0, 0, 1);

  seq.reset();
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid3.setIndex(idx);

    grid3.getSurfCoordSys(0, norm, tan1, tan2);
    Lucee::Vec3<double> lx(tan1), ly(tan2), lz(norm), t1t2;
    
    LC_ASSERT("Testing normal", ux.dot(lz) == 1.0);
    t1t2 = lx.cross(ly);
    for (unsigned i=0; i<3; ++i)
      LC_ASSERT("Testing right-handedness of coord-sys", t1t2[i] == lz[i]);
  }

  seq.reset();
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid3.setIndex(idx);

    grid3.getSurfCoordSys(1, norm, tan1, tan2);
    Lucee::Vec3<double> lx(tan1), ly(tan2), lz(norm), t1t2;
    
    LC_ASSERT("Testing normal", uy.dot(lz) == 1.0);
    t1t2 = lx.cross(ly);
    for (unsigned i=0; i<3; ++i)
      LC_ASSERT("Testing right-handedness of coord-sys", t1t2[i] == lz[i]);
  }

  seq.reset();
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid3.setIndex(idx);

    grid3.getSurfCoordSys(2, norm, tan1, tan2);
    Lucee::Vec3<double> lx(tan1), ly(tan2), lz(norm), t1t2;
    
    LC_ASSERT("Testing normal", uz.dot(lz) == 1.0);
    t1t2 = lx.cross(ly);
    for (unsigned i=0; i<3; ++i)
      LC_ASSERT("Testing right-handedness of coord-sys", t1t2[i] == lz[i]);
  }

}

int
main(int argc, char **argv)
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  LC_MPI_BEGIN_TESTS("lcrectcartgrid");
#else
  LC_BEGIN_TESTS("lcrectcartgrid");
#endif
  test_1();
  test_2();
  test_3();
#ifdef HAVE_MPI
  LC_MPI_END_TESTS;
#else
  LC_END_TESTS;
#endif
}
