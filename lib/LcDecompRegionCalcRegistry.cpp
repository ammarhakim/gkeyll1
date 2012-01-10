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

namespace Lucee
{
  void registerDecompRegionCalc(Lucee::LuaState& L)
  {
// register region decomp calculators
    new Lucee::ObjRegistry<Lucee::DecompRegionCalcIfc<1>, Lucee::CartProdDecompRegionCalc<1> >;
    new Lucee::ObjRegistry<Lucee::DecompRegionCalcIfc<2>, Lucee::CartProdDecompRegionCalc<2> >;
    new Lucee::ObjRegistry<Lucee::DecompRegionCalcIfc<3>, Lucee::CartProdDecompRegionCalc<3> >;
    new Lucee::ObjRegistry<Lucee::DecompRegionCalcIfc<4>, Lucee::CartProdDecompRegionCalc<4> >;
    new Lucee::ObjRegistry<Lucee::DecompRegionCalcIfc<5>, Lucee::CartProdDecompRegionCalc<5> >;
    new Lucee::ObjRegistry<Lucee::DecompRegionCalcIfc<6>, Lucee::CartProdDecompRegionCalc<6> >;
    new Lucee::ObjRegistry<Lucee::DecompRegionCalcIfc<7>, Lucee::CartProdDecompRegionCalc<7> >;
 } 
}

