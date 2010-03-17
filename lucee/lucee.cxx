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
#include <LcLogStream.h>
#include <LcLogger.h>
#include <LcLogger.h>
#include <LcLuaState.h>
#include <LcObjCreator.h>
#include <LcRegisterModules.h>
#include <LcSolverIfc.h>
#include <LcStreamHandler.h>

// std includes
#include <fstream>
#include <iostream>
#include <string>
#include <ctime>

int
main(int argc, char **argv)
{
  Lucee::CmdLineArgs cmdParser("lucee");
// add command line options
  cmdParser.addArg("i", "INPUT", "Input file");
  cmdParser.addArg("o", "OUTPUT-PREFIX", "Prefix for all output files");
  cmdParser.addArg("verbosity", "VERBOSITY", "Verbosity of log messages."
    " Should be one of disable,\n   debug, info, error. Defaults to info.");
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
  infoStrm << "Reading input file " << inpFile << std::endl;
  if (luaL_loadfile(L, inpFile.c_str()) || lua_pcall(L, 0, 0, 0))
  {
    std::cerr << "Error parsing input file " << inpFile << std::endl;
    exit(1);
  }

  double t0, t1;
// get simulation start time
  lua_getglobal(L, "tStart");
  if (! lua_type(L, -1) == LUA_TNUMBER)
    t0 = 0.0;
  else
    t0 = lua_tonumber(L, -1);
// get simulation end time
  lua_getglobal(L, "tEnd");
  if (! lua_type(L, -1) == LUA_TNUMBER)
    t1 = 1.0;
  else
    t1 = lua_tonumber(L, -1);

  if (t1 < t0)
  {
    std::cerr << "tEnd must be greater than tStart" << std::endl;
    exit(1);
  }

// get number of output frames to write
  unsigned outFrames;
  lua_getglobal(L, "outFrames");
  if (! lua_type(L, -1) == LUA_TNUMBER)
    outFrames = 1;
  else
    outFrames = lua_tonumber(L, -1);

// put top-level simulation table on stack
  lua_getglobal(L, "simulation");
// construct table object
  Lucee::LuaTable tbl(L, "simulation");

// get kind of simulation object
  std::string kind = tbl.getKind();
  infoStrm << "Creating top level simulation object '" << kind << "'" << std::endl;
// create a new simulation of this kind
  Lucee::SolverIfc *sim = Lucee::ObjCreator<Lucee::SolverIfc>
    ::getNew(kind);

  try
  {
// go through simulation boot-strap process
    sim->readInput(tbl);
    sim->buildData();
    sim->buildAlgorithms();

// check if we are restarting
    if (cmdParser.hasSwitch("r"))
      sim->restoreFromFile(inpFile);
    else
      sim->initialize();
    
// set current time
    sim->setCurrTime(t0);

// dump initial data before running simulation
    sim->writeToFile(outPrefix, 0);

    time_t start = time(0);
    struct tm *timeinfo;
    timeinfo = localtime(&start);
    infoStrm << std::endl;
    infoStrm << "Simulation started at time " << asctime(timeinfo) << std::endl;

// run simulation
    sim->advance(t1);

    time_t end = time(0); // time at end of main loop
    timeinfo = localtime ( &end );
    infoStrm << "Simulation finished at time " << asctime(timeinfo) << std::endl;

// dump final data after running simulation
    sim->writeToFile(outPrefix, outFrames);

// shut-down simulaiton
    sim->finalize();
  }
  catch (Lucee::Except& lce)
  {
    infoStrm << "Lucee exception ..." << std::endl;
    infoStrm << lce.what() << std::endl;
  }

  delete sim;

  return 0;
}
