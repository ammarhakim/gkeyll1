/**
 * @file	lcparticlecont.cxx
 *
 * @brief	Benchmarking particle containers.
 */

// lucee includes
#include <LcArray.h>
#include <LcConstFieldPtr.h>
#include <LcFieldPtr.h>
#include <LcField.h>
#include <LcParticleBase.h>
#include <LcParticleLayout.h>
#include <LcRegion.h>
#include <LcRowMajorIndexer.h>
#include <LcRowMajorSequencer.h>
#include <LcTest.h>

// boost includes
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>

// std includes
#include <ctime>
#include <iostream>
#include <list>
#include <vector>

// 3D/3V, double precision particle
class StandardParticle : public Lucee::ParticleBase<double, 
  Lucee::ParticleLayout<0, 3, 6> >
{
};

// 3D/3V, double precision particle, with extra marker flag to
// indicate if particle was moved from another cell
class StandardMarkedParticle : public Lucee::ParticleBase<double, 
  Lucee::ParticleLayout<0, 3, 6> >
{
  public:
    bool isNative() const { return isNativeFlag; }
    void setIsNative(bool inFlag) { isNativeFlag = inFlag; }
  private:
    bool isNativeFlag;
};

template <class PTCL>
class BorisPusher
{
  public:
    BorisPusher(double q, double m)
      : q(q), m(m)
    {
    }

    void push(double dt, double E[], double B[], PTCL& p)
    {
      double qmdt = 0.5*q/m*dt;
      double vm[3], vp[3], vprime[3], t[3], s[3], cp[3];

// compute t and s vectors
      for (unsigned d=0; d<3; ++d)
        t[d] = qmdt*B[d];
      double tNorm = t[0]*t[0] + t[1]*t[1] + t[2]*t[2];
      for (unsigned d=0; d<3; ++d)
        s[d] = 2*t[d]/(1+tNorm);
      
// half-step electric field update
      for (unsigned d=0; d<3; ++d)
        vm[d] = p.v(d) + qmdt*E[d];

// rotation around magnetic field
      crossProd(vm, t, cp);
      for (unsigned d=0; d<3; ++d)
        vprime[d] = vm[d] + cp[d];
      
      crossProd(vprime, s, cp);
      for (unsigned d=0; d<3; ++d)
        vp[d] = vm[d] + cp[d];
      
// half-step electric field update: this gives final particle velocity
      for (unsigned d=0; d<3; ++d)
        p.v(d) = vp[d]-qmdt*E[d];
      
// update particle position
      for (unsigned d=0; d<3; ++d)
        p.x(d) += dt*p.v(d);
    }

  private:
    void crossProd(const double a[], const double b[], double axb[])
    {
      axb[0] =  a[1]*b[2]-a[2]*b[1];
      axb[1] =  a[2]*b[0]-a[0]*b[2];
      axb[2] =  a[0]*b[1]-a[1]*b[0];
    }

    double q, m; // charge and mass
};

// particles stored as a global, unsorted linked list.
void
pushParticlesListNaive(
  std::list<StandardParticle>& ptclList,
  const Lucee::Region<2, double> domain,
  const Lucee::Field<2, double>& EM,
  unsigned nsteps, double dt)
{
// assumptions: EM fields are stored at nodes. Hence, the number of
// cells in each direction is one less than the number of nodes. No
// manual optimizations are done: i.e the implementation is maximally
// naive, something one would do as a first-cut

  BorisPusher<StandardParticle> pusher(-1.0, 1.0); // electrons

  Lucee::ConstFieldPtr<double> emPtr_11 = EM.createConstPtr();
  Lucee::ConstFieldPtr<double> emPtr_21 = EM.createConstPtr();
  Lucee::ConstFieldPtr<double> emPtr_12 = EM.createConstPtr();
  Lucee::ConstFieldPtr<double> emPtr_22 = EM.createConstPtr();

  Lucee::Region<2, int> rgn = EM.getRegion();

  double lower[2], upper[2], dx[2];
  for (unsigned d=0; d<2; ++d)
  {
    lower[d] = domain.getLower(d);
    upper[d] = domain.getUpper(d);
    dx[d] = (upper[d]-lower[d])/(rgn.getShape(d)-1); // numCells = numNodes-1
  }

  double EB[6], ploc[2];
  int idx[2];

  for (unsigned ns=0; ns<nsteps; ++ns)
  {
    //std::cout << ptclList.begin()->x(0) << " " << ptclList.begin()->x(1) << std::endl;

    std::list<StandardParticle>::iterator pItr = ptclList.begin();
    for ( ; pItr != ptclList.end(); ++pItr)
    {
// compute particle cell index
      for (unsigned d=0; d<2; ++d)
      {
        idx[d] = (int) (pItr->x(d)-lower[d])/dx[d];
        ploc[d] = (pItr->x(d) - (lower[d]+dx[d]*idx[d]) )/dx[d]; // incremental normalized position in cell
      }

// set pointers to EM field
      EM.setPtr(emPtr_11, idx[0], idx[1]);
      EM.setPtr(emPtr_21, idx[0]+1, idx[1]);
      EM.setPtr(emPtr_12, idx[0], idx[1]+1);
      EM.setPtr(emPtr_22, idx[0]+1, idx[1]+1);

// interpolate fields to particle location
      for (unsigned i=0; i<6; ++i)
      {
        double a1 = (1-ploc[0])*emPtr_11[i] + ploc[0]*emPtr_21[i];
        double a2 = (1-ploc[0])*emPtr_12[i] + ploc[0]*emPtr_22[i];
        EB[i] = (1-ploc[1])*a1 + ploc[1]*a2;
      }
// move particle
      pusher.push(dt, &EB[0], &EB[3], *pItr);

// check if particle has left domain (periodic BCs)
      for (unsigned d=0; d<2; ++d)
      {
        if (pItr->x(d) > upper[d])
          pItr->x(d) = lower[d] + (pItr->x(d)-upper[d]);
        if (pItr->x(d) < lower[d])
          pItr->x(d) = upper[d] - (lower[d]-pItr->x(d));
      }
    }
  }
}

// profile the above method to see where time is spent

// Store cell index and only increments in cell. However, this needs
// updating this information. Advantage: potentially better
// performance due to pre-computed values.

// Store particles in a per-cell list. This involves copying particles
// between lists. Potential gain: cache performance is improved.

void
testSingleParticlesListNaive()
{
  std::list<StandardParticle> ptclList;
// add a single particle, just to test
  StandardParticle p;
  p.x(0) = 1.0;
  p.x(1) = 0.0;
  p.x(2) = 0.0;

  p.v(0) = 0.0;
  p.v(1) = 5.0;
  p.v(2) = 0.0;

  ptclList.push_back(p);

  double qbm = 1.0;
  double dt = qbm/4.0;
  double pi = 3.141592654;
  unsigned nsteps = (unsigned) pi/dt;
  
  double dLo[2] = {-5, -5};
  double dUp[2] = {5, 5};
  Lucee::Region<2, double> domain(dLo, dUp);

// set field
  int shape[2] = {41, 41};
  Lucee::Region<2, int> grid(shape);
  Lucee::Field<2, double> EM(grid, 6, 0.0);
  Lucee::FieldPtr<double> emPtr = EM.createPtr();
  for (int i=grid.getLower(0); i<grid.getUpper(0); ++i)
  {
    for (int j=grid.getLower(1); j<grid.getUpper(1); ++j)
    {
      EM.setPtr(emPtr, i, j);
// electric field
      emPtr[0] = 0.0;
      emPtr[1] = 0.0;
      emPtr[2] = 0.0;
// magnetic field
      emPtr[3] = 0.0;
      emPtr[4] = 0.0;
      emPtr[5] = 1.0;
    }
  }

// push particles
  pushParticlesListNaive(ptclList, domain, EM, nsteps, dt);
}

void
testParticlesListNaive()
{
  std::list<StandardParticle> ptclList;

  double qbm = 1.0;
  double dt = qbm/10.0;
  double pi = 3.141592654;
  unsigned nsteps = 5000;
  unsigned npart = 40; // particles per cell
  
  double dLo[2] = {-5, -5};
  double dUp[2] = {5, 5};
  Lucee::Region<2, double> domain(dLo, dUp);

// set field
  int shape[2] = {41, 41};
  Lucee::Region<2, int> grid(shape);
  Lucee::Field<2, double> EM(grid, 6, 0.0);
  Lucee::FieldPtr<double> emPtr = EM.createPtr();
  for (int i=grid.getLower(0); i<grid.getUpper(0); ++i)
  {
    for (int j=grid.getLower(1); j<grid.getUpper(1); ++j)
    {
      EM.setPtr(emPtr, i, j);
// electric field
      emPtr[0] = 0.0;
      emPtr[1] = 0.0;
      emPtr[2] = 0.0;
// magnetic field
      emPtr[3] = 0.0;
      emPtr[4] = 0.0;
      emPtr[5] = 1.0;
    }
  }

  for (int i=grid.getLower(0); i<grid.getUpper(0)-1; ++i)
  {
    for (int j=grid.getLower(1); j<grid.getUpper(1)-1; ++j)
    {
      for (unsigned p=0; p<npart; ++p)
      {
        StandardParticle p;
// for now, just load all particles at cell center
        p.x(0) = domain.getLower(0) + 0.5*i;
        p.x(1) = domain.getLower(1) + 0.5*j;
        p.x(2) = 0.0;

        p.v(0) = 0.0;
        p.v(1) = 1.0;
        p.v(2) = 0.0;

        ptclList.push_back(p);
      }
    }
  }

// push particles
  clock_t start = clock();
  pushParticlesListNaive(ptclList, domain, EM, nsteps, dt);
  clock_t end = clock();
  double tsec = (end-start)/ (double) CLOCKS_PER_SEC;
  size_t numPushes = (grid.getShape(0)-1)*(grid.getShape(1)-1)*nsteps*npart;
  std::cout <<  tsec << " for flat-list.  "
            << " With " << numPushes << " pushes, resulting in "
            << tsec/numPushes << " per particle" << std::endl;
}

// Particle list in each cell
struct CellList
{
    std::list<StandardMarkedParticle> ptcls;
};


// Returns number of particles not in the cell they are supposed to
// be. This number should be exactly zero.
unsigned
checkIfAllPtclsInCell(const std::vector<CellList>& ptclList, 
  const Lucee::Region<2, double>& domain, 
  const Lucee::Region<2, int>& grid)
{
  double lower[2], upper[2], dx[2];
  for (unsigned d=0; d<2; ++d)
  {
    lower[d] = domain.getLower(d);
    upper[d] = domain.getUpper(d);
    dx[d] = (upper[d]-lower[d])/(grid.getShape(d)-1); // numCells = numNodes-1
  }

  int idx[2], intShape[2];
  intShape[0] = grid.getShape(0)-1;
  intShape[1] = grid.getShape(1)-1;
  Lucee::Region<2, int> interior(intShape);
  Lucee::RowMajorIndexer<2> indexer(interior);

  unsigned nOutOfBoundPtcls = 0;
  for (int i=interior.getLower(0); i<interior.getUpper(0); ++i)
  {
    for (int j=interior.getLower(1); j<interior.getUpper(1); ++j)
    {
      idx[0] = i; idx[1] = j;
      int cellIdx = indexer.getIndex(i,j);
      std::list<StandardMarkedParticle>::const_iterator pItr = ptclList[cellIdx].ptcls.begin();
      for ( ; pItr != ptclList[cellIdx].ptcls.end(); ++pItr)
      {
        int newIdx[2];
// ensure particle is in current cell
        for (unsigned d=0; d<2; ++d)
          newIdx[d] = (int) (pItr->x(d)-lower[d])/dx[d];
        if ((newIdx[0] != idx[0]) || (newIdx[1] != idx[1]))
          ++nOutOfBoundPtcls;
      }
    }
  }
  return nOutOfBoundPtcls;
}

// particles stored as a global, unsorted linked list.
void
pushParticlesPerCellList(
  std::vector<CellList>& ptclList,
  const Lucee::Region<2, double>& domain,
  const Lucee::Field<2, double>& EM,
  unsigned nsteps, double dt)
{
// assumptions: EM fields are stored at nodes. Hence, the number of
// cells in each direction is one less than the number of nodes.

  BorisPusher<StandardMarkedParticle> pusher(-1.0, 1.0); // electrons

  Lucee::ConstFieldPtr<double> emPtr_11 = EM.createConstPtr();
  Lucee::ConstFieldPtr<double> emPtr_21 = EM.createConstPtr();
  Lucee::ConstFieldPtr<double> emPtr_12 = EM.createConstPtr();
  Lucee::ConstFieldPtr<double> emPtr_22 = EM.createConstPtr();

  Lucee::Region<2, int> rgn = EM.getRegion();

// for timing
  double ptclPushTm = 0, ptclMoveTm = 0;

  double lower[2], upper[2], dx[2];
  for (unsigned d=0; d<2; ++d)
  {
    lower[d] = domain.getLower(d);
    upper[d] = domain.getUpper(d);
    dx[d] = (upper[d]-lower[d])/(rgn.getShape(d)-1); // numCells = numNodes-1
  }

  double EB[6], ploc[2];
  int idx[2];
  int intShape[2];
  intShape[0] = rgn.getShape(0)-1;
  intShape[1] = rgn.getShape(1)-1;
  Lucee::Region<2, int> interior(intShape);
  Lucee::RowMajorIndexer<2> indexer(interior);

  for (unsigned ns=0; ns<nsteps; ++ns)
  {
// main loop is over cells
    for (int i=interior.getLower(0); i<interior.getUpper(0); ++i)
    {
      for (int j=interior.getLower(1); j<interior.getUpper(1); ++j)
      {
        idx[0] = i; idx[1] = j;
// set pointers to EM field
        EM.setPtr(emPtr_11, idx[0], idx[1]);
        EM.setPtr(emPtr_21, idx[0]+1, idx[1]);
        EM.setPtr(emPtr_12, idx[0], idx[1]+1);
        EM.setPtr(emPtr_22, idx[0]+1, idx[1]+1);

// coordinates of lower-left and upper-right corners in cell
        double lowcc[2], upcc[2];
        for (unsigned d=0; d<2; ++d)
        {
          lowcc[d] = lower[d] + idx[d]*dx[d];
          upcc[d] = lower[d] + (idx[d]+1)*dx[d];
        }

        int cellIdx = indexer.getIndex(i,j);

// push particles current cell
        std::list<StandardMarkedParticle>::iterator pItr = ptclList[cellIdx].ptcls.begin();
        for ( ; pItr != ptclList[cellIdx].ptcls.end(); ++pItr)
        {
          if (pItr->isNative())
          { // only move if the particle is "native" to this cell
// compute incremental normalized position of particle in cell
            for (unsigned d=0; d<2; ++d)
              ploc[d] = (pItr->x(d)-lowcc[d])/dx[d];
          
// interpolate field to particle location
            for (unsigned i=0; i<6; ++i)
            {
              double a1 = (1-ploc[0])*emPtr_11[i] + ploc[0]*emPtr_21[i];
              double a2 = (1-ploc[0])*emPtr_12[i] + ploc[0]*emPtr_22[i];
              EB[i] = (1-ploc[1])*a1 + ploc[1]*a2;
            }
// move particle
            pusher.push(dt, &EB[0], &EB[3], *pItr);

// check if particle has left domain (periodic BCs)
            for (unsigned d=0; d<2; ++d)
            {
              if (pItr->x(d) > upper[d])
                pItr->x(d) = lower[d] + (pItr->x(d)-upper[d]);
              if (pItr->x(d) < lower[d])
                pItr->x(d) = upper[d] - (lower[d]-pItr->x(d));
            }
          }
        }

// loop over particles again, and copy those that have moved to
// appropriate cell
        pItr = ptclList[cellIdx].ptcls.begin();
        for ( ; pItr != ptclList[cellIdx].ptcls.end(); )
        {
          if (pItr->isNative())
          {
            int newIdx[2];
// check if particle has moved out of current cell
            for (unsigned d=0; d<2; ++d)
              newIdx[d] = (int) (pItr->x(d)-lower[d])/dx[d];
            
            if ((newIdx[0] != idx[0]) || (newIdx[1] != idx[1]))
            { // particle moved to another cell, remove it
// add it to new location, marking it as "non-native"
              pItr->setIsNative(false);
              unsigned newCellIdx = indexer.getIndex(newIdx);
              ptclList[newCellIdx].ptcls.push_back(*pItr);
// erase it from here
              pItr = ptclList[cellIdx].ptcls.erase(pItr); // sets it to next element
            }
            else
              ++pItr; // bump
          }
          else
          { // make non-native particle, native
            pItr->setIsNative(true);
            ++pItr; // bump
          }
        }
      }
    }

// final loop, making all particles local
    for (int i=interior.getLower(0); i<interior.getUpper(0); ++i)
    {
      for (int j=interior.getLower(1); j<interior.getUpper(1); ++j)
      {
        int cellIdx = indexer.getIndex(i,j);
        std::list<StandardMarkedParticle>::iterator pItr = ptclList[cellIdx].ptcls.begin();
        for ( ; pItr != ptclList[cellIdx].ptcls.end(); ++pItr)
          pItr->setIsNative(true);
      }
    }    
  }
  std::cout << "Particle push time " << ptclPushTm
            << ". Particle move time " << ptclMoveTm << std::endl;
}

void
testParticlesPerCellList()
{
  double qbm = 1.0;
  double dt = qbm/10.0;
  double pi = 3.141592654;
  unsigned nsteps = 5000;
  unsigned npart = 40; // particles per cell
  
  double dLo[2] = {-5, -5};
  double dUp[2] = {5, 5};
  Lucee::Region<2, double> domain(dLo, dUp);

// set field
  int shape[2] = {41, 41};
  Lucee::Region<2, int> grid(shape);
  Lucee::Field<2, double> EM(grid, 6, 0.0);
  Lucee::FieldPtr<double> emPtr = EM.createPtr();
  for (int i=grid.getLower(0); i<grid.getUpper(0); ++i)
  {
    for (int j=grid.getLower(1); j<grid.getUpper(1); ++j)
    {
      EM.setPtr(emPtr, i, j);
// electric field
      emPtr[0] = 0.0;
      emPtr[1] = 0.0;
      emPtr[2] = 0.0;
// magnetic field
      emPtr[3] = 0.0;
      emPtr[4] = 0.0;
      emPtr[5] = 1.0;
    }
  }

// allocate memory for particles
  unsigned ncells = (grid.getShape(0)-1)*(grid.getShape(1)-1);
  std::vector<CellList> ptclList(ncells);

  int intShape[2];
  intShape[0] = grid.getShape(0)-1;
  intShape[1] = grid.getShape(1)-1;
  Lucee::Region<2, int> interior(intShape);
  Lucee::RowMajorIndexer<2> indexer(interior);

  for (int i=grid.getLower(0); i<grid.getUpper(0)-1; ++i)
  {
    for (int j=grid.getLower(1); j<grid.getUpper(1)-1; ++j)
    {
      int idx = indexer.getIndex(i,j); // cell index
      for (unsigned p=0; p<npart; ++p)
      {
        StandardMarkedParticle p;
// for now, just load all particles at cell center
        p.x(0) = domain.getLower(0) + 0.5*i;
        p.x(1) = domain.getLower(1) + 0.5*j;
        p.x(2) = 0.0;

        p.v(0) = 0.0;
        p.v(1) = 1.0;
        p.v(2) = 0.0;

        p.setIsNative(true); // initially all particles are "native"

        ptclList[idx].ptcls.push_back(p);
      }
    }
  }

// push particles
  clock_t start = clock();
  pushParticlesPerCellList(ptclList, domain, EM, nsteps, dt);
  clock_t end = clock();
  double tsec = (end-start)/ (double) CLOCKS_PER_SEC;
  size_t numPushes = (grid.getShape(0)-1)*(grid.getShape(1)-1)*nsteps*npart;
  std::cout <<  tsec << " for local-list. "
            << " With " << numPushes << " pushes, resulting in "
            << tsec/numPushes << " per particle" << std::endl;
  unsigned noob = checkIfAllPtclsInCell(ptclList, domain, grid);
  std::cout << "Number of out-of-bound particles (should be zero) " 
            << noob << std::endl;
}

int
main(int argc, char *argv[])
{
  //testParticlesListNaive();
  testParticlesPerCellList();
  return 0;
}

