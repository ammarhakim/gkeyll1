/**
 * @file	gkeyll.cxx
 *
 * @brief	Top level driver for gkeyll simulations
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

// boost includes
#include <boost/lexical_cast.hpp>

// txbase includes
#ifdef HAVE_MPI
# include <TxMpiBase.h>
#ifdef HAVE_FFTW3
# include <fftw3-mpi.h>
#endif
#else
# include <TxSelfBase.h>
#endif

//  petsc includes
#ifdef HAVE_PETSC
#include <petsc.h>
#endif

// std includes
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

int
main(int argc, char **argv)
{
// initialize MPI if building in parallel
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#ifdef HAVE_FFTW3
  fftw_mpi_init();
#endif
#endif

  TxCommBase *comm = Loki::SingletonHolder<Lucee::Globals>
    ::Instance().comm;

  Lucee::CmdLineArgs cmdParser("lucee");
  cmdParser.addArg("i", "INPUT", "Input file");
  cmdParser.addArg("o", "OUTPUT-PREFIX", "Prefix for all output files");
  cmdParser.addArg("p", "PETSC-OPTIONS", "Options database file for PetSc");
  cmdParser.addArg("verbosity", "VERBOSITY", "Verbosity of log messages."
    " Should be one of disabled,\n   debug, info, error. Defaults to info.");
  cmdParser.addArg("baseVerbosity", "VERBOSITY", "Verbosity of log messages from each processor."
    " Should be one of disabled,\n   debug, info, error. Defaults to debug.");
  cmdParser.addArg("r", "RESTART-FRAME", "Frame to restart from.");

  cmdParser.parse(argc, argv);
  if (cmdParser.hasSwitch("h")) 
  {
    cmdParser.showHelp();
    exit(1);
  }

  std::string inpFile;
  if (cmdParser.hasArg("i"))
    inpFile = cmdParser.getArg("i");
  else
  {
    if (comm->getRank() == 0)
    {
      std::cerr << "** Input file not specified" << std::endl;
      cmdParser.showHelp();
    }
    exit(1);
  }

  bool isRestarting = false;
  int rFrame = 0;
  if (cmdParser.hasArg("r"))
  {
    std::string reStr = cmdParser.getArg("r");
    isRestarting = true;
    reStr = cmdParser.getArg("r");
    rFrame = boost::lexical_cast<int>(reStr);
  }

#ifdef HAVE_PETSC
// initialize PetSc
  std::string ptscFile;
  if (cmdParser.hasArg("p"))
  { // database file specified
    ptscFile = cmdParser.getArg("p");
    std::ifstream ptscFl(ptscFile.c_str());
    if (! ptscFl )
    {
      std::cerr << "Unable to open PetSc options file " << ptscFile << std::endl;
      exit(1);
    }
    PetscInitialize(&argc, &argv, ptscFile.c_str(), PETSC_NULL);
  }
  else
  {
    PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
  }
#endif

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
  if (cmdParser.hasArg("baseVerbosity"))
    logger.setLevel(cmdParser.getArg("baseVerbosity"));
  else
    logger.setLevel("debug");
// create file stream
  std::ostringstream snm;
  snm << outPrefix << "_" << comm->getRank() << ".log";
  Lucee::FileHandler fhndlr(snm.str());
  fhndlr.attachToLogger("lucee");

// create console logger
  Lucee::Logger& conLogger = Lucee::Logger::create("lucee.console");
  if (cmdParser.hasArg("verbosity"))
    conLogger.setLevel(cmdParser.getArg("verbosity"));
  else
    conLogger.setLevel("info");
// create console stream
  Lucee::StreamHandler conStrm(std::cout);
  if (comm->getRank() == 0)
// write messages to console only on rank 0
    conStrm.attachToLogger("lucee.console");

// get stream for logging
  Lucee::LogStream infoStrm = conLogger.getInfoStream();

// load input file using Lua
  Lucee::LuaState L;
// set global state
  Loki::SingletonHolder<Lucee::Globals>
    ::Instance().L = &L;

// load lua library: this must be done before loading input file
  Lucee::registerModules(L);

  bool failed = false; // flag to indicate if run failed

// add command line options to the top-level module
  static const luaL_Reg topFuncs[] = { {NULL, NULL} };
  luaL_register(L, "Lucee", topFuncs);

  lua_pushstring(L, inpFile.c_str());
  lua_setfield(L, -2, "InputFile");

  lua_pushstring(L, outPrefix.c_str());
  lua_setfield(L, -2, "OutPrefix");

  lua_pushboolean(L, isRestarting);
  lua_setfield(L, -2, "IsRestarting");

  lua_pushnumber(L, rFrame);
  lua_setfield(L, -2, "RestartFrame");

  lua_pop(L, 1); // done adding command line stuff

  infoStrm << "** This is Gkeyll!" << std::endl;
  time_t start = time(0); // time at start of main loop
  clock_t start_t = clock();
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
    failed = true;
  }

  time_t end = time(0); // time at end of main loop
  clock_t end_t = clock();
  timeinfo = localtime ( &end );
  infoStrm << "Simulation took " << (double) (end_t-start_t)/CLOCKS_PER_SEC
           << " seconds and finished at time " << asctime(timeinfo) << std::endl;

#ifdef HAVE_PETSC
  PetscFinalize();
#endif

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  exit(0); // PROBABLY NOT A GOOD THING TO DO
  return failed ? 1 : 0;
}
