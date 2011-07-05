/**
 * @file	lucee.cxx
 *
 * @brief	Top level driver for all Lucee simulations
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// lucee includes
#include <LcCmdLineArgs.h>
#include <LcFileHandler.h>
#include <LcGlobals.h>
#include <LcLogStream.h>
#include <LcLogger.h>
#include <LcLuaState.h>
#include <LcRegisterModules.h>
#include <LcSolverIfc.h>
#include <LcStreamHandler.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

int
main(int argc, char **argv)
{
  Lucee::CmdLineArgs cmdParser("lucee");
// add command line options
  cmdParser.addArg("i", "INPUT", "Input file");
  cmdParser.addArg("o", "OUTPUT-PREFIX", "Prefix for all output files");
  cmdParser.addArg("verbosity", "VERBOSITY", "Verbosity of log messages."
    " Should be one of disabled,\n   debug, info, error. Defaults to info.");
  cmdParser.addSwitch("r", "Restart simulation");

// parse command line
  cmdParser.parse(argc, argv);
// show help if requested
  if (cmdParser.hasSwitch("h")) 
  {
    cmdParser.showHelp();
    exit(1);
  }

// determine input file
  std::string inpFile;
  if (cmdParser.hasArg("i"))
    inpFile = cmdParser.getArg("i");
  else
  {
    std::cerr << "** Input file not specified" << std::endl;
    cmdParser.showHelp();
    exit(1);
  }

// check if input file exist
  std::ifstream inp(inpFile.c_str());
  if (! inp )
  {
    std::cerr << "Unable to open file " << inpFile << std::endl;
    exit(1);
  }

// create output prefix
  std::string outPrefix;
  if (cmdParser.hasArg("o"))
    outPrefix = cmdParser.getArg("o");
  else
  {
// use input file name sans the .lua extension
    std::string snm = inpFile;
    unsigned trunc = inpFile.find_last_of(".", snm.size());
    if (trunc > 0)
      snm.erase(trunc, snm.size());
    outPrefix = snm;
  }
// store output prefix in globals
  Loki::SingletonHolder<Lucee::Globals>
    ::Instance().outPrefix = outPrefix;

// create top-level logger
  Lucee::Logger& logger = Lucee::Logger::create("lucee");
  logger.setLevel("debug"); // base logger should log everything
// create file stream
  Lucee::FileHandler fhndlr(outPrefix + "_0.log");
  fhndlr.attachToLogger("lucee");

// create console logger
  Lucee::Logger& conLogger = Lucee::Logger::create("lucee.console");
  if (cmdParser.hasArg("verbosity"))
    conLogger.setLevel(cmdParser.getArg("verbosity"));
  else
    conLogger.setLevel("info");
// create console stream
  Lucee::StreamHandler conStrm(std::cout);
  conStrm.attachToLogger("lucee.console");

// get stream for logging
  Lucee::LogStream infoStrm = conLogger.getInfoStream();

// load input file using Lua
  Lucee::LuaState L;
// load lua library: this must be done before loading input file
  Lucee::registerModules(L);

  infoStrm << "** Welcome to Lucee!" << std::endl;

  time_t start = time(0); // time at start of main loop
  struct tm * timeinfo;
  timeinfo = localtime ( &start );
  infoStrm << "Simulation started at time " << asctime(timeinfo) << std::endl;
  try 
  {
    if (luaL_loadfile(L, inpFile.c_str()) || lua_pcall(L, 0, 0, 0))
    {
      std::cerr << "Error parsing input file: " << inpFile << std::endl;
      std::string err(lua_tostring(L, -1));
      lua_pop(L, 1);
      std::cerr << err << std::endl;
      exit(1);
    }
  }
  catch (Lucee::Except& lce)
  {
    infoStrm << "*** Lucee exception!" << std::endl;
    infoStrm << lce.what() << std::endl;
  }

  time_t end = time(0); // time at end of main loop
  timeinfo = localtime ( &end );
  infoStrm << "Simulation finished at time " << asctime(timeinfo) << std::endl;

  return 0;
}
