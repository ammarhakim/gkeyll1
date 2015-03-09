/**
 * @file	LcBasicObj.cpp
 *
 * @brief	Interface class for basic Lucee objects.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBasicObj.h>

// std includes
#include <string>

namespace Lucee
{
/** Check if communicator is valid */
  static bool _isValid(TxCommBase *c) { return (bool) c; }

  BasicObj::BasicObj()
    : nm("--NO-NAME--") 
  {
  }

  BasicObj::BasicObj(const std::string& nm)
    : nm(nm) 
  {
  }

  BasicObj::~BasicObj() 
  {
  }

  void
  BasicObj::readInput(Lucee::LuaTable& tbl)
  {
  }

  void
  BasicObj::initialize()
  {
  }

  void
  BasicObj::setComm(TxCommBase* ci)
  {
    comm = ci;
  }

  TxCommBase*
  BasicObj::getComm() const
  {
    return comm;
  }

  bool
  BasicObj::isValidOnRank() const
  {
// THIS NEEDS TO CHANGE
    return true; 
  }

  void
  BasicObj::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
  }

  std::string
  BasicObj::getName() const { return nm; }

  void
  BasicObj::setName(const std::string& name) 
  {
    nm = name;
  }
}
