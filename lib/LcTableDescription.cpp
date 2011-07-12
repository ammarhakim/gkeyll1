/**
 * @file	LcTableDescription.cpp
 *
 * @brief	Description of a single value in a Lua table.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcTableDescription.h>

// std includes
#include <sstream>

namespace Lucee
{
  TableDescription::TableDescription(const std::string& nm)
    : name(nm)
  {
  }

  void
  TableDescription::checkAndSet(Lucee::LuaTable& tbl)
  {
    bool pass = true;
    std::ostringstream errMsg;
// check and set values
    pass &= checkAndSetValues<int>(tbl, errMsg);
    pass &= checkAndSetValues<double>(tbl, errMsg);
    pass &= checkAndSetValues<std::string>(tbl, errMsg);
// check and set vectors
    pass &= checkAndSetVectors<int>(tbl, errMsg);
    pass &= checkAndSetVectors<double>(tbl, errMsg);
    pass &= checkAndSetVectors<std::string>(tbl, errMsg);

    if (!pass)
      throw Lucee::Except(errMsg.str());
  }
}
