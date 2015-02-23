/**
 * @file	lcparsum.cxx
 *
 * @brief	Basic Map-Reduce algorithm
 */

// MPI includes
#include <mpi.h>

// std includes
#include <cstdlib>
#include <iostream>
#include <vector>


int
main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

// size of array
  unsigned nsize = 10000;

// number of PUs
  int npu;
  MPI_Comm_size(MPI_COMM_WORLD, &npu);
// rank of current PU
  int rnk;
  MPI_Comm_rank(comm, &rnk);

// store start index and size of each sub-domain
  std::vector<unsigned> domStartIdx(npu);
  std::vector<unsigned> domShape(npu);

// decompose global domain into npu sub-domains
  int numVals = nsize/npu;
  int extraVals = nsize % npu;

// initially, each domain is the same size
  for (unsigned i=0; i<npu; ++i)
    domShape[i] = numVals;
// distribute remaining cells so domains are almost balanced
  for (unsigned i=0; i<extraVals; ++i)
    domShape[i] += 1;
// set initial index for each sub-domain
  domStartIdx[0] = 0;
  for (unsigned i=1; i<npu; ++i)
    domStartIdx[i] = domStartIdx[i-1]+domShape[i-1];

// allocate space for field
  std::vector<double> field(domShape[rnk]);

// initialize field
  for (unsigned i=0; i<domShape[rnk]; ++i)
    field[i] = 1.0;

// compute local sum of squares
  double sqSum = 0.0;
  for (unsigned i=0; i<domShape[rnk]; ++i)
    sqSum += field[i]*field[i];

// now, sum across domains
  double totalSum;
  MPI_Allreduce(&sqSum, &totalSum, 1, MPI_DOUBLE, MPI_SUM, comm);

// print out result only from rank 0
  if (rnk == 0)
    std::cout << "Sum is " << totalSum << std::endl;

  MPI_Finalize();
  return 0;
}
