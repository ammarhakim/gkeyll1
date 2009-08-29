/**
 * @file	LcSettableObject.cc
 *
 * @brief	Base class for all Lucee objects which can be set from name/values pairs.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// lucee includes
#include <LcExcept.h>
#include <LcSettableObject.h>

namespace Lucee
{
  SettableObject::SettableObject(const std::string& name)
    : name(name) 
  {
  }

  SettableObject::~SettableObject()
  {
  }

  std::string
  SettableObject::getName() const 
  {
    return name;
  }

  void
  SettableObject::addData(const std::string& nm, int *data, const std::string& help)
  {
    SettableData<int> dat;
    dat.data = data;
    dat.help = help;
    intMap[nm] = dat;
  }

  void
  SettableObject::addData(const std::string& nm, double *data, const std::string& help)
  {
    SettableData<double> dat;
    dat.data = data;
    dat.help = help;
    doubleMap[nm] = dat;
  }

  void
  SettableObject::addData(const std::string& nm, std::string *data, const std::string& help)
  {
    SettableData<std::string> dat;
    dat.data = data;
    dat.help = help;
    strMap[nm] = dat;
  }

  void
  SettableObject::setData(const std::string& nm, int val)
  {
    std::map<std::string, SettableData<int> >::iterator itr
      = intMap.find(nm);
    if (itr == intMap.end())
    {
      Lucee::Except lex("SettableObject::setData: Int data ");
      lex << nm << " does not exist in object " << name << std::endl;
      throw lex;
    }
    *(itr->second.data) = val;
  }

  void
  SettableObject::setData(const std::string& nm, double val)
  {
    std::map<std::string, SettableData<double> >::iterator itr
      = doubleMap.find(nm);
    if (itr == doubleMap.end())
    {
      Lucee::Except lex("SettableObject::setData: Double data ");
      lex << nm << " does not exist in object " << name << std::endl;
      throw lex;
    }
    *(itr->second.data) = val;
  }

  void
  SettableObject::setData(const std::string& nm, const std::string& val)
  {
    std::map<std::string, SettableData<std::string> >::iterator itr
      = strMap.find(nm);
    if (itr == strMap.end())
    {
      Lucee::Except lex("SettableObject::setData: String data ");
      lex << nm << " does not exist in object " << name << std::endl;
      throw lex;
    }
    *(itr->second.data) = val;
  }
}
