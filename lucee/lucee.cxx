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
#include <LcLuaState.h>
#include <LcSimulation.h>

// std includes
#include <fstream>
#include <iostream>
#include <string>

int
main(int argc, char **argv)
{
  Lucee::CmdLineArgs cmdParser("lucee");
// add command line options
  cmdParser.addArg("i", "INPUT", "Input file");
  cmdParser.addArg("op", "OUTPUT-PREFIX", "Prefix for all output files");
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

// create top-level simulation object
  Lucee::Simulation sim;

// load lua library: this must be done before loading input file

// load input file using Lua to determine some global variables
  Lucee::LuaState L;
  if (luaL_loadfile(L, inpFile.c_str()) || lua_pcall(L, 0, 0, 0))
  {
    std::cerr << "Error parsing input file " << inpFile << std::endl;
    exit(1);
  }

  double t0, t1;

// get start and end times
  lua_getglobal(L, "tStart");
  if (! lua_type(L, -1) == LUA_TNUMBER)
    t0 = 0.0;
  else
    t0 = lua_tonumber(L, -1);

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

// put the top-level simulation table on stack
  lua_getglobal(L, "simulation");
// construct table object with this table
  Lucee::LuaTable tbl(L, "simulation");

// go through simulation boot-strap process
  sim.readInput(tbl);
  sim.buildData();
  sim.buildAlgorithms();

// check if we are restarting
  if (cmdParser.hasSwitch("r"))
    sim.restoreFromFile(inpFile);
  else
    sim.initialize();

// set current time
  sim.setCurrTime(t0);
// run simulation
  sim.advance(t1);

// shut-down simulaiton
  sim.finalize();

  return 0;
}
