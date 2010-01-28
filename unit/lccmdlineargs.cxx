/**
 * @file	lccmdlineargs.cxx
 *
 * @brief	Unit tests parsing command lines.
 *
 * @version	$Id: lccmdlineargs.cxx 162 2009-08-28 20:07:10Z a.hakim777 $
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// lucee includes
#include <LcTest.h>
#include <LcCmdLineArgs.h>

int
main(void)
{
  LC_BEGIN_TESTS("lccmdlineargs");

// create command line object
  Lucee::CmdLineArgs cmd("lucee");
  cmd.addSwitch("r", "Restart simulation");
  cmd.addSwitch("dumph5", "Dump HDF5");
  cmd.addArg("i", "INPUT", "Input file");
  cmd.addArg("o", "OUTPUT", "Output prefix");
  cmd.addArg("email", "EMAIL", "Send results to this email address");

// print help message
  cmd.showHelp();

  LC_END_TESTS;
}
