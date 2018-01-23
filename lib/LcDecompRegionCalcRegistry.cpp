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
#include <LcCartGeneralDecompRegionCalc.h>
#include <LcCartProdDecompRegionCalc.h>
#include <LcDecompRegionCalcRegistry.h>
#include <LcRegisteredObjList.h>
#include <LcSubCartProdDecompRegionCalc.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  void registerDecompRegionCalc(Lucee::LuaState& L)
  {
// register region decomp calculators
    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::DecompRegionCalcIfc<1> > >
      ::Instance()
      .append<Lucee::CartProdDecompRegionCalc<1> >()
      .append<Lucee::CartGeneralDecompRegionCalc<1> >()

      .append<Lucee::SubCartProdDecompRegionCalc<1,2> >()
      .append<Lucee::SubCartProdDecompRegionCalc<1,3> >()
      .append<Lucee::SubCartProdDecompRegionCalc<1,4> >()
      .append<Lucee::SubCartProdDecompRegionCalc<1,5> >()
      .append<Lucee::SubCartProdDecompRegionCalc<1,6> >();
      

    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::DecompRegionCalcIfc<2> > >
      ::Instance()
      .append<Lucee::CartProdDecompRegionCalc<2> >()
      .append<Lucee::CartGeneralDecompRegionCalc<2> >()

      .append<Lucee::SubCartProdDecompRegionCalc<2,3> >()
      .append<Lucee::SubCartProdDecompRegionCalc<2,4> >()
      .append<Lucee::SubCartProdDecompRegionCalc<2,5> >()
      .append<Lucee::SubCartProdDecompRegionCalc<2,6> >();

    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::DecompRegionCalcIfc<3> > >
      ::Instance()
      .append<Lucee::CartProdDecompRegionCalc<3> >()
      .append<Lucee::CartGeneralDecompRegionCalc<3> >()

      .append<Lucee::SubCartProdDecompRegionCalc<3,4> >()
      .append<Lucee::SubCartProdDecompRegionCalc<3,5> >()
      .append<Lucee::SubCartProdDecompRegionCalc<3,6> >();

    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::DecompRegionCalcIfc<4> > >
      ::Instance()
      .append<Lucee::CartProdDecompRegionCalc<4> >()
      .append<Lucee::CartGeneralDecompRegionCalc<4> >()

      .append<Lucee::SubCartProdDecompRegionCalc<4,5> >() 
      .append<Lucee::SubCartProdDecompRegionCalc<4,6> >();
    
    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::DecompRegionCalcIfc<5> > >
      ::Instance()
      .append<Lucee::CartProdDecompRegionCalc<5> >()
      .append<Lucee::CartGeneralDecompRegionCalc<5> >();

    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::DecompRegionCalcIfc<6> > >
      ::Instance()
      .append<Lucee::CartProdDecompRegionCalc<6> >()
      .append<Lucee::CartGeneralDecompRegionCalc<6> >();
  }
}

