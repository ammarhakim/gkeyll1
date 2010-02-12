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
#include <LcCmdLineArgs.h>
#include <LcExcept.h>
#include <LcTest.h>

int
main(void)
{
  LC_BEGIN_TESTS("lccmdlineargs");

// create command line object
  Lucee::CmdLineArgs cmd("lucee");
  cmd.addSwitch("r", "Restart simulation");
  cmd.addSwitch("xml", "Dump XML");
  cmd.addArg("i", "INPUT", "Input file");
  cmd.addArg("o", "OUTPUT", "Output prefix");
  cmd.addArg("email", "EMAIL", "Send results to this email address");

// print help message
  //cmd.showHelp();

// test it
  char *argv1[] = 
    {
      "lccmdlineargs",
      "-r",
      "-xml",
      "-i",
      "input-file",
      "-o",
    "output-file"
    };
  cmd.parse(7, argv1);

  LC_ASSERT("Checking for switch", cmd.hasSwitch("r") == true);
  LC_ASSERT("Checking for switch", cmd.hasSwitch("xml") == true);
  LC_ASSERT("Checking for argument", cmd.hasArg("i") == true);
  LC_ASSERT("Checking for argument value", cmd.getArg("i") == "input-file");
  LC_ASSERT("Checking for argument value", cmd.getArg("o") == "output-file");

// test it
  char *argv2[] = 
    {
      "lccmdlineargs",
      "-xml",
      "-i",
      "input-file",
    };
  cmd.parse(4, argv2);

  LC_ASSERT("Checking for switch", cmd.hasSwitch("r") == false);
  LC_ASSERT("Checking for switch", cmd.hasSwitch("xml") == true);
  LC_ASSERT("Checking for argument", cmd.hasArg("i") == true);
  LC_ASSERT("Checking for argument value", cmd.getArg("i") == "input-file");
  LC_ASSERT("Checking for argument", cmd.hasArg("o") == false);
  LC_RAISES("Checking if exception is thrown", cmd.getArg("o"), Lucee::Except);

// test it
  char *argv3[] = 
    {
      "lccmdlineargs",
      "-x",
      "-y",
      "unknown",
      "-xml",
      "-i",
      "input-file"
    };
  cmd.parse(7, argv3);

  LC_ASSERT("Checking for switch", cmd.hasSwitch("r") == false);
  LC_ASSERT("Checking for switch", cmd.hasSwitch("xml") == true);
  LC_ASSERT("Checking for argument", cmd.hasArg("i") == true);
  LC_ASSERT("Checking for argument value", cmd.getArg("i") == "input-file");
  LC_ASSERT("Checking for argument", cmd.hasArg("o") == false);
  LC_RAISES("Checking if exception is thrown", cmd.getArg("o"), Lucee::Except);

// test it
  char *argv4[] = 
    {
      "lccmdlineargs",
      "-xml",
      "file-1",
      "file-2",
      "file-3"
    };
  cmd.parse(5, argv4);
  
// fetch extra arguments
  std::vector<std::string> extraArgs = cmd.getExtraArgs();
  
  LC_ASSERT("Checking if switch was specified", cmd.hasSwitch("xml") == true);
  LC_ASSERT("Checking if switch was specified", cmd.hasSwitch("r") == false);
  LC_ASSERT("Checking if number of extra arguments is correct",
    extraArgs.size() == 3);
  LC_ASSERT("Checking if extra argument is correct",
    extraArgs[0] == "file-1");
  LC_ASSERT("Checking if extra argument is correct",
    extraArgs[1] == "file-2");
  LC_ASSERT("Checking if extra argument is correct",
    extraArgs[2] == "file-3");

  LC_END_TESTS;
}
