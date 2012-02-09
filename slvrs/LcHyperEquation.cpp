/**
 * @file	LcHyperEquation.cpp
 *
 * @brief	Interface to hyperbolic equations.
 */
// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcExcept.h>
#include <LcHyperEquation.h>
#include <LcPointerHolder.h>


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
  HyperEquation::primitive(const Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double>& v) const
  {
    throw Lucee::Except("HyperEquation::primitive: Method not implemented");
  }

  void
  HyperEquation::conserved(const Lucee::ConstFieldPtr<double>& v, Lucee::FieldPtr<double>& q) const
  {
    throw Lucee::Except("HyperEquation::conserved: Method not implemented");
  }

  void
  HyperEquation::waves(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& jump, 
    const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
    Lucee::Matrix<double>& waves, Lucee::FieldPtr<double>& s)
  {
    throw Lucee::Except("HyperEquation::waves: Method not implemented");
  }

  double
  HyperEquation::numericalFlux(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
    Lucee::FieldPtr<double>& f)
  {
    throw Lucee::Except("HyperEquation::numericalFlux: Method not implemented");
    return 0; // should never come here
  }

  void
  HyperEquation::eigensystem(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& q,
    Lucee::Vector<double>& ev, Lucee::Matrix<double>& rev, Lucee::Matrix<double>& lev)
  {
    throw Lucee::Except("HyperEquation::eigensystem: Method not implemented");
  }

  void
  HyperEquation::quasiLinearMatrix(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& v, Lucee::Matrix<double>& qlMat)
  {
    throw Lucee::Except("HyperEquation::quasiLinearMatrix: Method not implemented");
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

  void
  HyperEquation::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
    lfm.appendFunc("primitive", luaPrimitive);
    lfm.appendFunc("conserved", luaConserved);
  }

  int
  HyperEquation::luaPrimitive(lua_State *L)
  {
    HyperEquation *hyp
      = Lucee::PointerHolder<HyperEquation>::getObjAsBase(L);

    if (lua_type(L, 2) != LUA_TUSERDATA)
    {
      Lucee::Except lce(
        "HyperEquation::luaPrimitive: Must provide a conserved variables field to 'primitive' method");
      throw lce;
    }
    Lucee::PointerHolder<Lucee::BasicObj> *cv =
      (Lucee::PointerHolder<Lucee::BasicObj>*) lua_touserdata(L, 2);

    if (lua_type(L, 3) != LUA_TUSERDATA)
    {
      Lucee::Except lce(
        "HyperEquation::luaPrimitive: Must provide a primitive variables field to 'primitive' method");
      throw lce;
    }
    Lucee::PointerHolder<Lucee::BasicObj> *pv =
      (Lucee::PointerHolder<Lucee::BasicObj>*) lua_touserdata(L, 3);

// call proper method based on type of object supplied
    if (cv->pointer->getType() == typeid(Lucee::StructGridField<1, double>).name()) {
      if (pv->pointer->getType() != cv->pointer->getType())
        throw Lucee::Except(
          "HyperEquation::luaPrimitive: Mismatched object types in 'primitive' method");
      hyp->calcPrimVars<1>(*(Lucee::StructGridField<1, double>*) cv->pointer, *(Lucee::StructGridField<1, double>*)pv->pointer);
    }
    else if (cv->pointer->getType() == typeid(Lucee::StructGridField<2, double>).name()) {
      if (pv->pointer->getType() != cv->pointer->getType())
        throw Lucee::Except(
          "HyperEquation::luaPrimitive: Mismatched object types in 'primitive' method");
      hyp->calcPrimVars<2>(*(Lucee::StructGridField<2, double>*) cv->pointer, *(Lucee::StructGridField<2, double>*)pv->pointer);
    }
    else if (cv->pointer->getType() == typeid(Lucee::StructGridField<3, double>).name()) {
      if (pv->pointer->getType() != cv->pointer->getType())
        throw Lucee::Except(
          "HyperEquation::luaPrimitive: Mismatched object types in 'primitive' method");
      hyp->calcPrimVars<3>(*(Lucee::StructGridField<3, double>*) cv->pointer, *(Lucee::StructGridField<3, double>*)pv->pointer);
    }
    else {
      throw Lucee::Except("HyperEquation::luaPrimitive: Incorrect object in method 'primitive'");
    }

    return 0;
  }

  int
  HyperEquation::luaConserved(lua_State *L)
  {
    HyperEquation *hyp
      = Lucee::PointerHolder<HyperEquation>::getObjAsBase(L);

    if (lua_type(L, 2) != LUA_TUSERDATA)
    {
      Lucee::Except lce(
        "HyperEquation::luaConserved: Must provide a conserved variables field to 'conserved' method");
      throw lce;
    }
    Lucee::PointerHolder<Lucee::BasicObj> *pv =
      (Lucee::PointerHolder<Lucee::BasicObj>*) lua_touserdata(L, 2);

    if (lua_type(L, 3) != LUA_TUSERDATA)
    {
      Lucee::Except lce(
        "HyperEquation::luaConserved: Must provide a conserved variables field to 'conserved' method");
      throw lce;
    }
    Lucee::PointerHolder<Lucee::BasicObj> *cv =
      (Lucee::PointerHolder<Lucee::BasicObj>*) lua_touserdata(L, 3);

// call proper method based on type of object supplied
    if (cv->pointer->getType() == typeid(Lucee::StructGridField<1, double>).name()) {
      if (pv->pointer->getType() != cv->pointer->getType())
        throw Lucee::Except(
          "HyperEquation::luaConserved: Mismatched object types in 'conserved' method");
      hyp->calcConsVars<1>(*(Lucee::StructGridField<1, double>*) pv->pointer, *(Lucee::StructGridField<1, double>*) cv->pointer);
    }
    else if (cv->pointer->getType() == typeid(Lucee::StructGridField<2, double>).name()) {
      if (pv->pointer->getType() != cv->pointer->getType())
        throw Lucee::Except(
          "HyperEquation::luaConserved: Mismatched object types in 'conserved' method");
      hyp->calcConsVars<2>(*(Lucee::StructGridField<2, double>*) pv->pointer, *(Lucee::StructGridField<2, double>*) cv->pointer);
    }
    else if (cv->pointer->getType() == typeid(Lucee::StructGridField<3, double>).name()) {
      if (pv->pointer->getType() != cv->pointer->getType())
        throw Lucee::Except(
          "HyperEquation::luaConserved: Mismatched object types in 'conserved' method");
      hyp->calcConsVars<3>(*(Lucee::StructGridField<3, double>*) pv->pointer, *(Lucee::StructGridField<3, double>*) cv->pointer);
    }
    else {
      throw Lucee::Except("HyperEquation::luaConserved: Incorrect object in method 'conserved'");
    }

    return 0;
  }
}
