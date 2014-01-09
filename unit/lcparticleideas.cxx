/**
 * @file	lcparticleideas.cxx
 *
 * @brief	Ideas to test for particles in Gkeyll
 */

// lucee includes
#include <LcParticleBase.h>
#include <LcParticleLayout.h>
#include <LcParticleList.h>
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
#include <iostream>
#include <vector>
#include <ctime>

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

void
test_ptcl_container()
{
  Lucee::ParticleList<StandardParticle<double> > plist(50);
  unsigned nalloc = 0;
// start adding particles
  StandardParticle<double> p;
  for (unsigned i=0; i<100000; ++i)
  {
    p.setx(0, i+0.5);
    p.setx(1, i+1.5);
    p.setx(2, i+2.5);

    p.setv(0, 2*i+0.5);
    p.setv(1, 2*i+1.5);
    p.setv(2, 2*i+2.5);

    if (plist.excessCapacity() <= 0)
      nalloc++;

    plist.addParticle(p);
  }

  unsigned i=127;
  LC_ASSERT("Testing particle", plist[i].x(0) == i+0.5);
  LC_ASSERT("Testing particle", plist[i].x(1) == i+1.5);
  LC_ASSERT("Testing particle", plist[i].x(2) == i+2.5);

  LC_ASSERT("Testing particle", plist[i].v(0) == 2*i+0.5);
  LC_ASSERT("Testing particle", plist[i].v(1) == 2*i+1.5);
  LC_ASSERT("Testing particle", plist[i].v(2) == 2*i+2.5);

  i=9867;
  LC_ASSERT("Testing particle", plist[i].x(0) == i+0.5);
  LC_ASSERT("Testing particle", plist[i].x(1) == i+1.5);
  LC_ASSERT("Testing particle", plist[i].x(2) == i+2.5);

  LC_ASSERT("Testing particle", plist[i].v(0) == 2*i+0.5);
  LC_ASSERT("Testing particle", plist[i].v(1) == 2*i+1.5);
  LC_ASSERT("Testing particle", plist[i].v(2) == 2*i+2.5);
}

struct IntParticle : public Lucee::ParticleBase<int,
  Lucee::ParticleLayout<0,0,1> >
{
    int getValue()
    { return this->getLoc(0); }
    void setValue(int val)
    { return this->setLoc(0, val); }
};

void
test_ptcl_update()
{
  typedef Lucee::ParticleList<StandardParticle<double> > PList_t;
  typedef PList_t::iterator PItr_t;
  typedef boost::minstd_rand base_generator_type;
  base_generator_type generator(42);

  int lower[2] = {0, 0};
  int upper[2] = {100, 100};
  Lucee::Region<2, int> dom(lower, upper);
  Lucee::RowMajorIndexer<2> domIdxr(dom);
  Lucee::RowMajorSequencer<2> seq(dom);

  double xlo = 0.0, ylo = 0.0;
  double xup = 1.0, yup = 1.0;
  double dx = (xup-xlo)/dom.getShape(0);
  double dy = (yup-ylo)/dom.getShape(1);
  double vmax = 1.0;

  boost::uniform_real<> uni_dist(0,1);
  boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);

  boost::uniform_real<> univ_dist(-vmax,vmax);
  boost::variate_generator<base_generator_type&, boost::uniform_real<> > 
    univ(generator, univ_dist);

  unsigned numPtclPerCell = 100;
// allocate space for particles
  std::vector<PList_t> particles(dom.getVolume());

  int idx[2];
// add particles
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    unsigned ci = domIdxr.getIndex(idx);

// bottom left node
    double xn = idx[0]*dx;
    double yn = idx[1]*dy;
    for (unsigned i=0; i<numPtclPerCell; ++i)
    {
      StandardParticle<double> ptcl;
      ptcl.setx(0, xn + uni());
      ptcl.setx(1, yn + uni());

      ptcl.setv(0, univ());
      ptcl.setv(1, univ());

      particles[ci].addParticle(ptcl);
    }
  }

  double cfl = 0.5;
  double dt = cfl*dx/vmax;
// update particles
  seq.reset();
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    unsigned ci = domIdxr.getIndex(idx);

    for (PItr_t pItr = particles[ci].begin(); pItr != particles[ci].end(); ++pItr)
    {
      pItr->incrx(0, dt*pItr->v(0));
      pItr->incrx(1, dt*pItr->v(1));
    }
  }
}

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

// rotation around magentic field
  crossProd(vm, t, cp);
  for (unsigned d=0; d<3; ++d)
    vprime[d] = vm[d] + cp[d];

  crossProd(vprime, s, cp);
  for (unsigned d=0; d<3; ++d)
    vp[d] = vm[d] + cp[d];

// half-step electric field update: this gives final particle velocity
  for (unsigned d=0; d<3; ++d)
    p.setv(d, vp[d]-qmdt*E[d]);

// update particle position
  for (unsigned d=0; d<3; ++d)
    p.incrx(d, dt*p.v(d));
}

int
main(int argc, char *argv[])
{
  LC_BEGIN_TESTS("lcparticleideas");
  test_ptcl_classes();
  test_ptcl_container();
  test_ptcl_update();

  LC_END_TESTS;
}

