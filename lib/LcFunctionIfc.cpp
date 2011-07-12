/**
 * @file	LcFunctionIfc.cpp
 *
 * @brief	Interface class for function objects.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcExcept.h>
#include <LcFunctionIfc.h>

namespace Lucee
{
  const char *FunctionIfc::id = "Function";

  FunctionIfc::FunctionIfc()
    : Lucee::BasicObj(FunctionIfc::id), inpSz(1), outSz(1)
  {
  }

  FunctionIfc::FunctionIfc(unsigned ni, unsigned no)
    : Lucee::BasicObj(FunctionIfc::id), inpSz(ni), outSz(no)
  {
  }

  void
  FunctionIfc::readInput(Lucee::LuaTable& tbl)
  {
    throw Lucee::Except("FunctionIfc::readInput: Method not implemented.");
  }
}
