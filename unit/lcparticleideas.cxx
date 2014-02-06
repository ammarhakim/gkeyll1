/**
 * @file	lcparticleideas.cxx
 *
 * @brief	Ideas to test for particles in Gkeyll
 */

// lucee includes
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
#include <boost/intrusive/list.hpp>

// std includes
#include <ctime>
#include <iostream>
#include <list>
#include <vector>

// 3D/3V
template <typename REAL>
class StandardParticle : public Lucee::ParticleBase<REAL, 
  Lucee::ParticleLayout<0, 3, 6> >
{
};

// 3D/3V
template <typename REAL>
struct WeightedParticle : public Lucee::ParticleBase<REAL, 
  Lucee::ParticleLayout<0, 3, 7> >
{
    REAL w() const
    { return this->getLoc(6); }

    void setw(REAL w)
    { this->setLoc(6, w); }
};

void
test_ptcl_classes()
{
  StandardParticle<float> p1;
  p1.setx(0, 1.0);
  p1.setx(1, 2.0);
  p1.setx(2, 3.0);

  p1.setv(0, 1.5);
  p1.setv(1, 2.5);
  p1.setv(2, 3.5);

  LC_ASSERT("Testing particle", p1.x(0) == 1.0);
  LC_ASSERT("Testing particle", p1.x(1) == 2.0);
  LC_ASSERT("Testing particle", p1.x(2) == 3.0);

  LC_ASSERT("Testing particle", p1.v(0) == 1.5);
  LC_ASSERT("Testing particle", p1.v(1) == 2.5);
  LC_ASSERT("Testing particle", p1.v(2) == 3.5);

  StandardParticle<float> p2;
  p2.setx(0, 2.0);
  p2.setx(1, 3.0);
  p2.setx(2, 4.0);

  p2.setv(0, 2.5);
  p2.setv(1, 3.5);
  p2.setv(2, 4.5);

  p1.copyParticle(p2);

  LC_ASSERT("Testing particle", p1.x(0) == 2.0);
  LC_ASSERT("Testing particle", p1.x(1) == 3.0);
  LC_ASSERT("Testing particle", p1.x(2) == 4.0);

  LC_ASSERT("Testing particle", p1.v(0) == 2.5);
  LC_ASSERT("Testing particle", p1.v(1) == 3.5);
  LC_ASSERT("Testing particle", p1.v(2) == 4.5);

  WeightedParticle<float> wp;
  wp.setx(0, 1.0);
  wp.setx(1, 2.0);
  wp.setx(2, 3.0);

  wp.setv(0, 1.5);
  wp.setv(1, 2.5);
  wp.setv(2, 3.5);

  wp.setw(5.6);

  LC_ASSERT("Testing particle", wp.x(0) == 1.0);
  LC_ASSERT("Testing particle", wp.x(1) == 2.0);
  LC_ASSERT("Testing particle", wp.x(2) == 3.0);

  LC_ASSERT("Testing particle", wp.v(0) == 1.5);
  LC_ASSERT("Testing particle", wp.v(1) == 2.5);
  LC_ASSERT("Testing particle", wp.v(2) == 3.5);

  LC_ASSERT("Testing particle", wp.w() == (float) 5.6);
}

struct IntParticle : public Lucee::ParticleBase<int,
  Lucee::ParticleLayout<0,0,1> >
{
    int getValue()
    { return this->getLoc(0); }
    void setValue(int val)
    { return this->setLoc(0, val); }
};

// The copying of particles in the following is to avoid memory
// overhead when deleting particles in the middle of a list. The
// assumption is that copying is cheap, while reallocations expensive.

// Another related implementation could be to use a linked-list
// (std::list) to store the particles. The extra copy would be
// eliminated, but then particles would no longer be contiguous in
// memory. This could be a problem with cache-misses.

void
deleteAndAddPtcl(std::vector<StandardParticle<double> >& ptclSrc, unsigned loc,
  std::vector<StandardParticle<double> >& ptclDest)
{
  ptclDest.push_back(ptclSrc[loc]);
  if (ptclSrc.size()>1)
    ptclSrc[loc] = ptclSrc[ptclSrc.size()-1];
  ptclSrc.pop_back();
}

void
deletePtcl(std::vector<StandardParticle<double> >& ptclSrc, unsigned loc)
{
  if (ptclSrc.size()>1)
    ptclSrc[loc] = ptclSrc[ptclSrc.size()-1];
  ptclSrc.pop_back();
}

void
test_ptcl_copy_delete()
{
  std::vector<StandardParticle<double> > ptclsCell1;
  std::vector<StandardParticle<double> > ptclsCell2;
  std::vector<StandardParticle<double> > ptclsCell3;

// push some particles in each list
  for (unsigned i=0; i<10; ++i)
  {
    StandardParticle<double> p;
    p.x(0) = i+0.5;
    p.v(0) = 2*i+0.5;
    ptclsCell1.push_back(p);
  }

  for (unsigned i=0; i<5; ++i)
  {
    StandardParticle<double> p;
    p.x(0) = i;
    p.v(0) = 2*i;
    ptclsCell2.push_back(p);
  }

  for (unsigned i=0; i<1; ++i)
  {
    StandardParticle<double> p;
    p.x(0) = i+0.25;
    p.v(0) = 2*i+0.25;
    ptclsCell3.push_back(p);
  }

  // std::cout << "In cell 1 (" << ptclsCell1.size() << ")" << std::endl;
  // for (unsigned i=0; i<ptclsCell1.size(); ++i)
  // {
  //   std::cout << " " << ptclsCell1[i].x(0) << " " << ptclsCell1[i].v(0) << std::endl;
  // }
  // std::cout << "In cell 2 ("  << ptclsCell2.size() << ")" << std::endl;
  // for (unsigned i=0; i<ptclsCell2.size(); ++i)
  // {
  //   std::cout << " " << ptclsCell2[i].x(0) << " " << ptclsCell2[i].v(0) << std::endl;
  // }

// now delete some stuff
  deleteAndAddPtcl(ptclsCell1, 5, ptclsCell2);
  deleteAndAddPtcl(ptclsCell3, 0, ptclsCell2);

// check it
  LC_ASSERT("Testing delete/add", ptclsCell1.size() == 9);
  LC_ASSERT("Testing delete/add", ptclsCell1[5].x(0) == 9+0.5);
  LC_ASSERT("Testing delete/add", ptclsCell1[5].v(0) == 2*9+0.5);

  LC_ASSERT("Testing delete/add", ptclsCell2.size() == 7);
  LC_ASSERT("Testing delete/add", ptclsCell2[5].x(0) == 5+0.5);
  LC_ASSERT("Testing delete/add", ptclsCell2[5].v(0) == 2*5+0.5);

  LC_ASSERT("Testing delete/add", ptclsCell2[6].x(0) == 0+0.25);
  LC_ASSERT("Testing delete/add", ptclsCell2[6].v(0) == 0*5+0.25);

  LC_ASSERT("Testing delete/add", ptclsCell3.size() == 0);

  // std::cout << "In cell 1 (" << ptclsCell1.size() << ")" << std::endl;
  // for (unsigned i=0; i<ptclsCell1.size(); ++i)
  // {
  //   std::cout << " " << ptclsCell1[i].x(0) << " " << ptclsCell1[i].v(0) << std::endl;
  // }
  // std::cout << "In cell 2 ("  << ptclsCell2.size() << ")" << std::endl;
  // for (unsigned i=0; i<ptclsCell2.size(); ++i)
  // {
  //   std::cout << " " << ptclsCell2[i].x(0) << " " << ptclsCell2[i].v(0) << std::endl;
  // }
  
}

// NOTE: For efficient particle-push, one needs to tell the compiler
// that pointers are not aliased. Or replace this with macro.

void crossProd(const double a[], const double b[], double axb[])
{
  axb[0] =  a[1]*b[2]-a[2]*b[1];
  axb[1] =  a[2]*b[0]-a[0]*b[2];
  axb[2] =  a[0]*b[1]-a[1]*b[0];
}

void
borisPush(double q, double m, double dt, StandardParticle<double>& p)
{
  double qmdt = 0.5*q/m*dt;
  double vm[3], vp[3], vprime[3], t[3], s[3], cp[3];
  double B[3] = {0.0, 0.0, 1.0};
  double E[3] = {0.0, 0.0, 0.0};

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

void
cyclotron()
{
  StandardParticle<double> p;
  p.setx(0, 1.0);
  p.setx(1, 0.0);
  p.setx(2, 0.0);

  p.setv(0, 0.0);
  p.setv(1, 1.0);
  p.setv(2, 0.0);

  std::cout << p.x(0) << " " << p.x(1) << std::endl;
  double qbm = 1.0;
  double dt = qbm/4.0;
  double pi = 3.141592654;
  unsigned nsteps = (unsigned) 2*pi/dt;
  for (unsigned i=0; i<3*nsteps; ++i)
  {
    borisPush(1.0, 1.0, dt, p);
    std::cout << p.x(0) << " " << p.x(1) << std::endl;
  }
}

void
cyclotronTime()
{
  unsigned npart = 10000;
  std::vector<StandardParticle<double> > p(npart);
  for (unsigned i=0; i<npart; ++i)
  {
    p[i].setx(0, std::rand());
    p[i].setx(1, 0.0);
    p[i].setx(2, 0.0);
    
    p[i].setv(0, 0.0);
    p[i].setv(1, std::rand());
    p[i].setv(2, 0.0);
  }

  unsigned ntries = 10000;
  double qbm = 1.0;
  double dt = qbm/4.0;
  double pi = 3.141592654;
  unsigned c = 0;
  for (unsigned t=0; t<ntries; ++t)
    for (unsigned pIdx=0; pIdx<npart; ++pIdx)
    {
      borisPush(1.0, 1.0, dt, p[pIdx]);
      c++;
    }
  std::cout << "Total pushes " << c << std::endl;
}

int
main(int argc, char *argv[])
{
  LC_BEGIN_TESTS("lcparticleideas");
  //test_ptcl_classes();
  //test_ptcl_copy_delete();
  //cyclotron();
  cyclotronTime();

  LC_END_TESTS;
}

