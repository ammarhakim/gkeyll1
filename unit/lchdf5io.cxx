/**
 * @file	lchdf5io.cxx
 *
 * @brief	Unit tests for Hdf5 txbase I/O classes.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcTest.h>

// txbase includes
#include <TxHdf5Base.h>
#include <TxIoBase.h>
#ifdef HAVE_MPI
# include <TxMpiBase.h>
#else
# include <TxSelfBase.h>
#endif

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
  TxCommBase* comm = new TxMpiBase();
#else
  TxCommBase* comm = new TxSelfBase();
#endif
  TxIoBase* io = new TxHdf5Base(comm);

  // open a file
  TxIoNodeType fn = io->createFile("hdf5io_test.h5");
  // write data to file
  TxIoNodeType dn = io->writeDataSet<double>(
    fn, "testdata", dataSetSize, dataSetBeg, dataSetLen, data);
  // write some attributes to this datanode
  io->writeAttribute(dn, "Time", 1.56);
  io->writeAttribute(dn, "PI", 3.141592654);

  // write a vector attribute
  std::vector<int> intVec;
  intVec.push_back(10);
  intVec.push_back(20);
  intVec.push_back(30);
  io->writeAttribute(dn, "Int Vect", intVec);

  std::vector<double> floatVec;
  floatVec.push_back(1.0);
  floatVec.push_back(2.0);
  floatVec.push_back(3.0);
  io->writeAttribute(dn, "Float Vect", floatVec);

  delete io;

  // now read it back in
  TxIoBase* ior = new TxHdf5Base(comm);
  TxIoNodeType fnr = ior->openFile("hdf5io_test.h5", "r");
  double *datar = new double[50];
  TxIoNodeType dnr = ior->readDataSet<double>(
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
  io->readAttribute(dnr, "Int Vect", intVecR);

  std::vector<double> floatVecR;
  io->readAttribute(dnr, "Float Vect", floatVecR);

  for (unsigned i=0; i<3; ++i)
  {
    LC_ASSERT("Checking for vector int attribute",
      intVecR[i]==intVec[i]);
    LC_ASSERT("Checking for vector float attribute",
      floatVecR[i]==floatVec[i]);
  }

  delete [] data;
  delete [] datar;
  delete ior;
}

int 
main(int argc, char *argv[])
{
  LC_BEGIN_TESTS("lchdf5io");
  test_1();
  LC_END_TESTS;
}
