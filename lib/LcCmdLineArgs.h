/**
 * @file	LcCmdLineArgs.h
 *
 * @brief	Class to parse command line arguments.
 *
 * @version	$Id: LcCmdLineArgs.h 222 2009-11-17 04:46:22Z a.hakim777 $
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

#ifndef LC_CMD_LINE_ARGS_H
#define LC_CMD_LINE_ARGS_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// std includes
#include <map>
#include <vector>
#include <string>

namespace Lucee
{
  class CmdLineArgs 
  {
    public:
/**
 * Create a new command line parse.
 *
 * @param name Name of program
 */
      CmdLineArgs(const std::string& name);

/**
 * Add a switch.
 *
 * @param nm Switch name.
 * @param help Help string.
 */
      void addSwitch(const std::string& nm, const std::string& help);

/**
 * Add an argument.
 *
 * @param nm Argument name.
 * @param disp Display name for help.
 * @param help Help string.
 */
      void addArg(const std::string& nm, const std::string& disp, const std::string& help);

/**
 * Show help message.
 */
      void showHelp() const;

/**
 * Parse command line argumanets specified.
 *
 * @param argc Number of arguments.
 * @param argv List of arguments.
 */
      void parse(int argc, char *argv[]);

/**
 * Check if a switch exists.
 *
 * @param nm Name of switch to check.
 * @return true if it exists, false otherwise.
 */
      bool hasSwitch(const std::string& nm) const;

/**
 * Check if an argument exists.
 *
 * @param nm Name of argument to check.
 * @return true if it exists, false otherwise.
 */
      bool hasArg(const std::string& nm) const;

/**
 * Return value for argument.
 *
 * @param nm Name of argument.
 * @return Value of argument.
 */
      std::string getArg(const std::string& nm) const;

/**
 * Return list of extra arguments.
 *
 * @return list of extra arguments.
 */
      std::vector<std::string> getExtraArgs() const 
      {
        return extraArgs;
      }

    private:
/**
 * Is this a switch?
 *
 * @param nm Name to test.
 * @return true if it is a switch, false otherwise.
 */
      bool isSwitch(const std::string& nm) const;

/**
 * Is this an argument?
 *
 * @param nm Name to test.
 * @return true if it is a name, false otherwise.
 */
      bool isArg(const std::string& nm) const;

/** Name of program */
      std::string progNm;
/** Map of switches to help string */
      std::map<std::string, std::string> switchMap;
/** Map of arguments to help string */
      std::map<std::string, std::string> argMap;
/** Map of arguments to display-string */
      std::map<std::string, std::string> dispMap;
/** Map of switches */
      std::map<std::string, bool> switchVal;
/** Map of arguments to values */
      std::map<std::string, std::string> nmVal;
/** List of extra arguments */
      std::vector<std::string> extraArgs;
  };
}

#endif // LC_CMD_LINE_ARGS_H


