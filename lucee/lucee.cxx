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
#include <string>
#include <iostream>

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

// create top-level simulation object
  Lucee::Simulation sim;

// now open input file to determine some global variables

// go through simulation boot-strap process
  sim.readInput();
  sim.buildData();
  sim.buildAlgorithms();

// check if we are restarting
  if (cmdParser.hasSwitch("r"))
    sim.restoreFromFile(inpFile);
  else
    sim.initialize();

// set current time
  sim.setCurrTime(0.0); // FROM INPUT FILE
// run simulation
  sim.advance(1.0); // FROM INPUT FILE

// shut-down simulaiton
  sim.finalize();

  return 0;
}
