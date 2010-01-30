/**
 * @file	LcCmdLineArgs.cpp
 *
 * @brief	Class to parse command line arguments.
 *
 * @version	$Id: LcCmdLineArgs.cpp 222 2009-11-17 04:46:22Z a.hakim777 $
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// std includes
#include <iostream>

// lucee includes
#include <LcCmdLineArgs.h>
#include <LcExcept.h>

namespace Lucee
{

  CmdLineArgs::CmdLineArgs(const std::string& name)
    : progNm(name) 
  {
    addSwitch("h", "Show help message and exit");
  }

  void
  CmdLineArgs::addSwitch(const std::string& nm, const std::string& help) 
  {
    switchMap[nm] = help;
  }

  void
  CmdLineArgs::addArg(const std::string& nm, const std::string& disp, const std::string& help) 
  {
    argMap[nm] = help;
    dispMap[nm] = disp;
  }

  void
  CmdLineArgs::showHelp() const 
  {
    std::cout << "Usage " << progNm << " OPTIONS extra-args" << std::endl;
    std::cout << " -h " << std::endl;
    std::cout << "   Show help message and exit" << std::endl;
    std::map<std::string, std::string>::const_iterator itr;
    for (itr=switchMap.begin(); itr!=switchMap.end(); ++itr) 
    {
// show help strings for all switches except for "-h"
      if (itr->first != "h") 
      {
        std::cout << " -" << itr->first << std::endl;
        std::cout << "  " << itr->second << std::endl;
      }
    }
    for (itr=argMap.begin(); itr!=argMap.end(); ++itr) 
    {
      std::map<std::string, std::string>::const_iterator dItr =
        dispMap.find(itr->first);
      std::cout << " -" << itr->first << " " <<  dItr->second << std::endl;
      std::cout << "   " << itr->second << std::endl;
    }
  }

  void
  CmdLineArgs::parse(int argc, char *argv[]) 
  {
// clear existing maps
    switchVal.clear();
    nmVal.clear();
    extraArgs.clear();

    int i=1;
// loop over each argument, putting it into proper maps
    while (i<argc) 
    {
      if (argv[i][0] == '-') 
      {
// this is an a switch or an argument
        char *nm = &argv[i][1]; // eat the '-'
        if (isSwitch(nm)) 
        {
          switchVal[nm] = true;
        }
        else if (isArg(nm)) 
        {
          char *val = argv[i+1];
          nmVal[nm] = val;
          i = i+1; // skip over value
        }
        i = i+1;
      }
      else 
      {
// this is an extra argument
        extraArgs.push_back(argv[i]);
        i = i+1;
      }
    }
  }

  bool
  CmdLineArgs::hasSwitch(const std::string& nm) const 
  {
    std::map<std::string, bool>::const_iterator itr
      = switchVal.find(nm);
    if (itr != switchVal.end())
      return true;
    return false;
  }

  bool
  CmdLineArgs::hasArg(const std::string& nm) const 
  {
    std::map<std::string, std::string>::const_iterator itr
      = nmVal.find(nm);
    if (itr != nmVal.end())
      return true;
    return false;
  }

  std::string
  CmdLineArgs::getArg(const std::string& nm) const 
  {
    std::map<std::string, std::string>::const_iterator itr
      = nmVal.find(nm);
    if (itr != nmVal.end())
      return itr->second;
    Lucee::Except lce("CmdLineArgs::getArg: Argument ");
    lce << nm << " does not exist" << std::endl;
    throw lce;
  }

  bool
  CmdLineArgs::isSwitch(const std::string& nm) const 
  {
    std::map<std::string, std::string>::const_iterator itr
      = switchMap.find(nm);
    if (itr != switchMap.end())
      return true;
    return false;
  }

  bool
  CmdLineArgs::isArg(const std::string& nm) const 
  {
    std::map<std::string, std::string>::const_iterator itr
      = argMap.find(nm);
    if (itr != argMap.end())
      return true;
    return false;
  }
}
