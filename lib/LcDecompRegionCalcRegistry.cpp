/**
 * @file LcDecompRegionCalcRegistry.cpp
 *
 * @brief	Method for registering basic library object.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCartProdDecompRegionCalc.h>
#include <LcDecompRegionCalcRegistry.h>
#include <LcRegisteredObjList.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  void registerDecompRegionCalc(Lucee::LuaState& L)
  {
// register region decomp calculators
    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::DecompRegionCalcIfc<1> > >
      ::Instance()
      .append<Lucee::CartProdDecompRegionCalc<1> >();

    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::DecompRegionCalcIfc<2> > >
      ::Instance()
      .append<Lucee::CartProdDecompRegionCalc<2> >();

    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::DecompRegionCalcIfc<3> > >
      ::Instance()
      .append<Lucee::CartProdDecompRegionCalc<3> >();

    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::DecompRegionCalcIfc<4> > >
      ::Instance()
      .append<Lucee::CartProdDecompRegionCalc<4> >();

    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::DecompRegionCalcIfc<5> > >
      ::Instance()
      .append<Lucee::CartProdDecompRegionCalc<5> >();
  }
}

