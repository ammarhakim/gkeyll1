/**
 * @file LcDecompRegionCalcRegistry.h
 *
 * @brief	Method for registering basic library object.
 */

#ifndef LC_DECOMP_REGION_CALC_REGISTRY_H
#define LC_DECOMP_REGION_CALC_REGISTRY_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcObjRegistry.h>

namespace Lucee
{
/**
 * Register objects and modules.
 */
  void registerDecompRegionCalc(Lucee::LuaState& L);
}

#endif // LC_DECOMP_REGION_CALC_REGISTRY_H
