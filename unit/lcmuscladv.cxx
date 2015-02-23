/**
 * @file	lcmuscladv.cxx
 *
 * @brief	2D advection equations using MUSCL agorithm
 */

// blitz inlcudes
#include <blitz/array.h>

// txbase includes
#include <mpi.h>

// std includes
#include <cstdlib>
#include <iostream>
#include <vector>

class Grid
{
  public:
/** Grid is defined by lower/upper bounds and number of cells */
    Grid(double lower[2], double upper[2], int cells[2]);

    double lower[2], upper[2];
    int cells[2];
};

Grid::Grid(double lower[2], double upper[2], int cells[2])
{
  for (unsigned d=0; d<2; ++d)
  {
    this->lower[d] = lower[d];
    this->upper[d] = upper[d];
    this->cells[d] = cells[d];
  }
}

class AdvectionUpdater
{
  public:
    AdvectionUpdater(const Grid& g, int cuts[2]);

/**
 * Advance solution by specified time-step
 */    
    void advance(double dt);

  private:
    Grid grid;
    int lowerIdx[2], upperIdx[2];

    void decompose(int cuts[2]);
};

AdvectionUpdater::AdvectionUpdater(const Grid& g, int cuts[2])
  : grid(g)
{
  this->decompose(cuts);
}

void AdvectionUpdater::decompose(int cuts[2])
{
}

int
main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);

  double lower[2] = {0.0, 0.0};
  double upper[2] = {1.0, 1.0};
  int cells[2] = {50, 50};
  Grid grid(lower, upper, cells);

  MPI_Finalize();
  return 0;
}
