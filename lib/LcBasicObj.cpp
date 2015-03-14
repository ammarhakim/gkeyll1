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
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <string>

namespace Lucee
{
/** Check if communicator is valid */
  static bool _isValid(TxCommBase *c) { return (bool) c; }

  BasicObj::BasicObj()
    : nm("--NO-NAME--") 
  {
// set communicators to global one
    this->setComm(Loki::SingletonHolder<Lucee::Globals>
      ::Instance().comm);
    setIsSafeToWrite(true);
  }

  BasicObj::BasicObj(const std::string& nm)
    : nm(nm) 
  {
// set communicators to global one
    this->setComm(Loki::SingletonHolder<Lucee::Globals>
      ::Instance().comm);
    setIsSafeToWrite(true);
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

  void
  BasicObj::setMomComm(TxCommBase* ci)
  {
    momComm = ci;
  }

  TxCommBase*
  BasicObj::getMomComm() const
  {
    return momComm;
  }

  void
  BasicObj::setIsSafeToWrite(bool stw)
  {
    safeToWrite = stw;
  }

  bool
  BasicObj::isSafeToWrite() const
  {
    return safeToWrite;
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
