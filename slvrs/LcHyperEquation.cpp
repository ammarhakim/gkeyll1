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
  HyperEquation::flux(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double>& f)
  {
    throw Lucee::Except("HyperEquation::flux: Method not implemented");
  }

  void
  HyperEquation::speeds(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& q, double s[2])
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
  HyperEquation::waves(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& jump, 
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
    eigensystem(c, qavg, ev, rev, lev);
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
  HyperEquation::eigensystem(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& q,
    Lucee::Vector<double>& ev, Lucee::Matrix<double>& rev, Lucee::Matrix<double>& lev)
  {
    throw Lucee::Except("HyperEquation::eigensystem: Method not implemented");
  }

  bool
  HyperEquation::isInvariantDomain(const Lucee::ConstFieldPtr<double>& q) const
  {
    throw Lucee::Except("HyperEquation::isInvariantDomain: Method not implemented");
  }

  void
  HyperEquation::qFluctuations(const Lucee::Matrix<double>& waves, const Lucee::FieldPtr<double>& s,
    Lucee::FieldPtr<double>& amdq, Lucee::FieldPtr<double>& apdq)
  {
    for (unsigned m=0; m<meqn; ++m)
    {
      amdq[m] = 0.0; apdq[m] = 0.0;
      for (unsigned mw=0; mw<mwave; ++mw)
      {
        if (s[mw] < 0.0)
        {
// left going wave
          amdq[m] += s[mw]*waves(m,mw);
        }
        else
        {
// right going waves
          apdq[m] += s[mw]*waves(m,mw);
        }
      }
    }
  }

  void
  HyperEquation::fFluctuations(const Lucee::Matrix<double>& waves, const Lucee::FieldPtr<double>& s,
    Lucee::FieldPtr<double>& amdq, Lucee::FieldPtr<double>& apdq)
  {
    for (unsigned m=0; m<meqn; ++m)
    {
      amdq[m] = 0.0; apdq[m] = 0.0;
      for (unsigned mw=0; mw<mwave; ++mw)
      {
        if (s[mw] < 0.0)
        {
// left going wave
          amdq[m] += waves(m,mw);
        }
        else if (s[mw] > 0.0)
        {
// right going waves
          apdq[m] += waves(m,mw);
        }
        else
        {
// zero wave: add contribution to apdq and amdq
          amdq[m] += 0.5*waves(m, mw);
          apdq[m] += 0.5*waves(m, mw);
        }
      }
    }
  }
}
