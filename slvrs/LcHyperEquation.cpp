/**
 * @file	LcHyperEquation.cpp
 *
 * @brief	Interface to hyperbolic equations.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */
// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcExcept.h>
#include <LcHyperEquation.h>

namespace Lucee
{
// set module name
  const char *HyperEquation::id = "HyperEquation";

  HyperEquation::HyperEquation(unsigned meqn, unsigned mwave)
    : meqn(meqn), mwave(mwave)
  {
  }

  void
  HyperEquation::readInput(Lucee::LuaTable& tbl)
  {
  }

  void
  HyperEquation::flux(const Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double>& f)
  {
    throw Lucee::Except("HyperEquation::flux: Method not implemented");
  }

  void
  HyperEquation::speeds(const Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double>& s)
  {
    throw Lucee::Except("HyperEquation::speeds: Method not implemented");
  }

  void
  HyperEquation::primitive(const Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double>& v)
  {
    throw Lucee::Except("HyperEquation::primitive: Method not implemented");
  }

  void
  HyperEquation::conserved(const Lucee::ConstFieldPtr<double>& v, Lucee::FieldPtr<double>& q)
  {
    throw Lucee::Except("HyperEquation::conserved: Method not implemented");
  }

  void
  HyperEquation::waves(const Lucee::ConstFieldPtr<double>& jump, 
    const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
    Lucee::Matrix<double>& waves, Lucee::FieldPtr<double>& s)
  {
// calculate left and right primitive variables
    Lucee::FieldPtr<double> vl(meqn), vr(meqn);
    primitive(ql, vl);
    primitive(qr, vr);
// average it
    Lucee::FieldPtr<double> vavg(meqn);
    for (unsigned i=0; i<meqn; ++i)
      vavg[i] = 0.5*(vl[i]+vr[i]);
// compute conserved variables
    Lucee::FieldPtr<double> qavg(meqn);
    conserved(vavg, qavg);

// compute eigensystem
    Lucee::Vector<double> ev(meqn);
    Lucee::Matrix<double> rev(meqn, meqn), lev(meqn, meqn);
    eigensystem(qavg, ev, rev, lev);
// split jump using left eigenvectors
    Lucee::Vector<double> alpha(meqn);
    for (unsigned p=0; p<meqn; ++p)
    {
      alpha[p] = 0.0;
      for (unsigned i=0; i<meqn; ++i)
        alpha[p] += lev(i,p)*jump[i];
    }
// compute waves
    for (unsigned p=0; p<meqn; ++p)
    {
      for (unsigned i=0; i<meqn; ++i)
        waves(i,p) = alpha[p]*rev(i,p);
      s[p] = ev[p]; // set speeds
    }
  }

  void
  HyperEquation::eigensystem(const Lucee::ConstFieldPtr<double>& q,
    Lucee::Vector<double>& ev, Lucee::Matrix<double>& rev, Lucee::Matrix<double>& lev)
  {
    throw Lucee::Except("HyperEquation::eigensystem: Method not implemented");
  }

  bool
  HyperEquation::isInvariantDomain(const Lucee::ConstFieldPtr<double>& q) const
  {
    throw Lucee::Except("HyperEquation::isInvariantDomain: Method not implemented");
  }
}
