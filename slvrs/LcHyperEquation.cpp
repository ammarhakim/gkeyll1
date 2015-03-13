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
#include <LcGlobals.h>
#include <LcHyperEquation.h>
#include <LcPointerHolder.h>

// loki includes
#include <loki/Singleton.h>

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
  HyperEquation::rotateToLocal(const Lucee::RectCoordSys& c, const double* inQ, double* outQ)
  {
    throw Lucee::Except("HyperEquation::rotateToLocal: Method not implemented");
  }

  void
  HyperEquation::rotateToGlobal(const Lucee::RectCoordSys& c, const double* inQ, double* outQ)
  {
    throw Lucee::Except("HyperEquation::rotateToGlobal: Method not implemented");
  }

  void
  HyperEquation::flux(const Lucee::RectCoordSys& c, const double* q, 
    const std::vector<const double*>& auxVars, double* f)
  {
    throw Lucee::Except("HyperEquation::flux: Method not implemented");
  }

  void
  HyperEquation::speeds(const Lucee::RectCoordSys& c, const double* q, double s[2])
  {
    throw Lucee::Except("HyperEquation::speeds: Method not implemented");
  }

  void
  HyperEquation::primitive(const double* q, double* v) const
  {
    throw Lucee::Except("HyperEquation::primitive: Method not implemented");
  }

  void
  HyperEquation::conserved(const double* v, double* q) const
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
    const double* ql, const double* qr, 
    const std::vector<const double*>& auxVarsl, const std::vector<const double*>& auxVarsr,
    double* f)
  {
    throw Lucee::Except("HyperEquation::numericalFlux: Method not implemented");
    return 0; // should never come here
  }

  void
  HyperEquation::projectOnLeftEigenvectors(const Lucee::RectCoordSys& c,
    const double* q, const double* vec, double* coeff)
  {
    throw Lucee::Except("HyperEquation::projectOnLeftEigenvectors: Method not implemented");
  }

  void
  HyperEquation::reconWithRightEigenvectors(const Lucee::RectCoordSys& c,
    const double* q, const double* coeff, double* vec)
  {
    throw Lucee::Except("HyperEquation::reconWithRightEigenvectors: Method not implemented");
  }

  void
  HyperEquation::eigensystem(const Lucee::RectCoordSys& c,
    const double *q,
    Lucee::Vector<double>& ev, Lucee::Matrix<double>& rev, Lucee::Matrix<double>& lev)
  {
    throw Lucee::Except("HyperEquation::eigensystem: Method not implemented");
  }

  void
  HyperEquation::quasiLinearMatrix(const Lucee::RectCoordSys& c,
    const double *v, Lucee::Matrix<double>& qlMat)
  {
    throw Lucee::Except("HyperEquation::quasiLinearMatrix: Method not implemented");
  }

  bool
  HyperEquation::isInvariantDomain(const double* q) const
  {
    throw Lucee::Except("HyperEquation::isInvariantDomain: Method not implemented");
  }

  void
  HyperEquation::qFluctuations(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
    const Lucee::Matrix<double>& waves, const Lucee::FieldPtr<double>& s,
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
  HyperEquation::fFluctuations(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
    const Lucee::Matrix<double>& waves, const Lucee::FieldPtr<double>& s,
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
    lfm.appendFunc("checkInvariantDomain", luaIsInvariantDomain);
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
    else if (cv->pointer->getType() == typeid(Lucee::StructGridField<2, double>).name()) 
    {
      if (pv->pointer->getType() != cv->pointer->getType())
        throw Lucee::Except(
          "HyperEquation::luaPrimitive: Mismatched object types in 'primitive' method");
      hyp->calcPrimVars<2>(*(Lucee::StructGridField<2, double>*) cv->pointer, *(Lucee::StructGridField<2, double>*)pv->pointer);
    }
    else if (cv->pointer->getType() == typeid(Lucee::StructGridField<3, double>).name()) 
    {
      if (pv->pointer->getType() != cv->pointer->getType())
        throw Lucee::Except(
          "HyperEquation::luaPrimitive: Mismatched object types in 'primitive' method");
      hyp->calcPrimVars<3>(*(Lucee::StructGridField<3, double>*) cv->pointer, *(Lucee::StructGridField<3, double>*)pv->pointer);
    }
    else 
    {
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
    if (cv->pointer->getType() == typeid(Lucee::StructGridField<1, double>).name()) 
    {
      if (pv->pointer->getType() != cv->pointer->getType())
        throw Lucee::Except(
          "HyperEquation::luaConserved: Mismatched object types in 'conserved' method");
      hyp->calcConsVars<1>(*(Lucee::StructGridField<1, double>*) pv->pointer, *(Lucee::StructGridField<1, double>*) cv->pointer);
    }
    else if (cv->pointer->getType() == typeid(Lucee::StructGridField<2, double>).name()) 
    {
      if (pv->pointer->getType() != cv->pointer->getType())
        throw Lucee::Except(
          "HyperEquation::luaConserved: Mismatched object types in 'conserved' method");
      hyp->calcConsVars<2>(*(Lucee::StructGridField<2, double>*) pv->pointer, *(Lucee::StructGridField<2, double>*) cv->pointer);
    }
    else if (cv->pointer->getType() == typeid(Lucee::StructGridField<3, double>).name()) 
    {
      if (pv->pointer->getType() != cv->pointer->getType())
        throw Lucee::Except(
          "HyperEquation::luaConserved: Mismatched object types in 'conserved' method");
      hyp->calcConsVars<3>(*(Lucee::StructGridField<3, double>*) pv->pointer, *(Lucee::StructGridField<3, double>*) cv->pointer);
    }
    else 
    {
      throw Lucee::Except("HyperEquation::luaConserved: Incorrect object in method 'conserved'");
    }

    return 0;
  }

  int
  HyperEquation::luaIsInvariantDomain(lua_State *L)
  {
    HyperEquation *hyp
      = Lucee::PointerHolder<HyperEquation>::getObjAsBase(L);

    if (lua_type(L, 2) != LUA_TUSERDATA)
    {
      Lucee::Except lce(
        "HyperEquation::luaIsInvariantDomain: Must provide a conserved variables field to 'checkInvariantDomain' method");
      throw lce;
    }
    Lucee::PointerHolder<Lucee::BasicObj> *cv =
      (Lucee::PointerHolder<Lucee::BasicObj>*) lua_touserdata(L, 2);

    int result = true;
// call proper method based on type of object supplied
    if (cv->pointer->getType() == typeid(Lucee::StructGridField<1, double>).name()) 
    {
      result = hyp->checkInvariantDomain<1>(*(Lucee::StructGridField<1, double>*) cv->pointer);
    }
    else if (cv->pointer->getType() == typeid(Lucee::StructGridField<2, double>).name()) 
    {
      result = hyp->checkInvariantDomain<2>(*(Lucee::StructGridField<2, double>*) cv->pointer);
    }
    else if (cv->pointer->getType() == typeid(Lucee::StructGridField<3, double>).name()) 
    {
      result = hyp->checkInvariantDomain<3>(*(Lucee::StructGridField<3, double>*) cv->pointer);
    }
    else 
    {
      throw Lucee::Except("HyperEquation::luaConserved: Incorrect object in method 'conserved'");
    }
// make sure all processors return same answer
    TxCommBase *comm = cv->pointer->getComm();
    int globalResult;
    comm->allreduce(1, &result, &globalResult, TX_AND);

// push results on stack
    if (globalResult)
      lua_pushboolean(L, 1);
    else
      lua_pushboolean(L, 0);

    return 1;
  }
}
