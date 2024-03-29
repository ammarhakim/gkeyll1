/**
 * @file	LcRtePhaseFunction.cpp
 *
 * @brief	Radiative transfer equation in homogeneous slab.
 */

// lucee includes
#include <LcPointerHolder.h>
#include <LcRtePhaseFunction.h>

namespace Lucee
{
// set class ID for use in registration system
  const char *RtePhaseFunction::id = "RtePhaseFunction";

  RtePhaseFunction::~RtePhaseFunction()
  {
  }

  void
  RtePhaseFunction::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
    lfm.appendFunc("print", luaPrintExpCoeffs);
  }

  int
  RtePhaseFunction::luaPrintExpCoeffs(lua_State *L)
  {
    RtePhaseFunction *pf 
      = Lucee::PointerHolder<RtePhaseFunction>::getObjAsBase(L);
    int nL = (int) lua_tonumber(L, 2);
    Lucee::Vector<double> coeffs = pf->getExpCoeffs(nL);
    for (unsigned i=0; i<coeffs.size(); ++i)
      std::cout << coeffs[i] << " ";
    std::cout << std::endl;

    return 0;
  }
}
