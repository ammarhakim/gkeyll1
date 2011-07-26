/**
 * @file	lchdf5io.cxx
 *
 * @brief	Unit tests for Lucee::Hdf5Io class
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcHdf5Io.h>
#include <LcTest.h>

// HDF5 includes
#include <hdf5.h>

void
test_1() 
{
  std::vector<size_t> dataSetSize;
  std::vector<size_t> dataSetBeg;
  std::vector<size_t> dataSetLen;
  int lower[] = {0, 0};
  int upper[] = {10, 5};

  for (unsigned i=0; i<2; ++i)
  {
    dataSetSize.push_back( upper[i] );
    dataSetBeg.push_back( lower[i] );
    dataSetLen.push_back( upper[i] );
  }

  // write data
  double *data = new double[50];
  unsigned i=0;
  for (int j=lower[0]; j<upper[0]; ++j)
    for (int k=lower[1]; k<upper[1]; ++k)
      data[i++] = j+10.0*k;

  // create io object
#ifdef HAVE_MPI
  Lucee::IoBase* io = new Lucee::Hdf5Io(MPI_COMM_WORLD, MPI_INFO_NULL);
#else
  Lucee::IoBase* io = new Lucee::Hdf5Io(0, 0);
#endif
  // open a file
  Lucee::IoNodeType fn = io->createFile("hdf5io_test.h5");
  // write data to file
  Lucee::IoNodeType dn = io->writeDataSet<double>(
    fn, "testdata", dataSetSize, dataSetBeg, dataSetLen, data);
  // write some attributes to this datanode
  io->writeAttribute(dn, "Time", 1.56);
  io->writeAttribute(dn, "PI", 3.141592654);

  // write a vector attribute
  std::vector<int> intVec;
  intVec.push_back(10);
  intVec.push_back(20);
  intVec.push_back(30);
  io->writeVecAttribute(dn, "Int Vect", intVec);

  std::vector<double> floatVec;
  floatVec.push_back(1.0);
  floatVec.push_back(2.0);
  floatVec.push_back(3.0);
  io->writeVecAttribute(dn, "Float Vect", floatVec);

  delete io;

  // now read it back in
#ifdef HAVE_MPI
  Lucee::IoBase* ior = new Lucee::Hdf5Io(MPI_COMM_WORLD, MPI_INFO_NULL);
#else
  Lucee::IoBase* ior = new Lucee::Hdf5Io(0, 0);
#endif
  Lucee::IoNodeType fnr = ior->openFile("hdf5io_test.h5", "r");
  double *datar = new double[50];
  Lucee::IoNodeType dnr = ior->readDataSet<double>(
    fnr, "testdata", dataSetBeg, dataSetLen, datar);
  // check if we read it in properly
  for (int j=0; j<50; ++j)
    LC_ASSERT("Testing if data read correctly",
      datar[j] == data[j]);

  double dtime, dpi;
  ior->readAttribute<double>(dnr, "Time", dtime);
  ior->readAttribute<double>(dnr, "PI", dpi);
  LC_ASSERT("Testing if time read properly",
    dtime == 1.56);
  LC_ASSERT("Testing if PI read properly",
    dpi ==3.141592654);

  // read a vector attribute
  std::vector<int> intVecR;
  //io->readVecAttribute(dnr, "Int Vect", intVecR);

  std::vector<double> floatVecR;
  //io->readVecAttribute(dnr, "Float Vect", floatVecR);

//   for (unsigned i=0; i<3; ++i)
//   {
//     LC_ASSERT("Checking for vector int attribute",
//       intVecR[i]==intVec[i]);
//     LC_ASSERT("Checking for vector float attribute",
//       floatVecR[i]==floatVec[i]);
//   }

  delete [] data;
  delete [] datar;
  delete ior;
}

int 
main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  LC_MPI_BEGIN_TESTS("lchdf5io");
#else
  LC_BEGIN_TESTS("lchdf5io");
#endif
  test_1();

#ifdef HAVE_MPI
  LC_MPI_END_TESTS;
  MPI_Finalize();
#else
  LC_END_TESTS;
#endif
  return 0;
}
