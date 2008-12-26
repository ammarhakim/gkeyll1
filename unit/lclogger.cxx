/**
 * @file	lclogger.cxx
 *
 * @brief	Unit tests for logging classes
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#include <lctest.h>
#include <lclogger.h>
#include <lclogstream.h>
#include <lcfilehandler.h>
#include <lcstreamhandler.h>

void
test_a()
{
  Lucee::Logger& l = Lucee::Logger::get("lucee");

  // set level for logger
  l.setLevel("debug");
  LC_ASSERT("Testing level", l.getLevelStr() == "debug");

  // create a stream handler
  std::ostringstream myStrm;
  Lucee::StreamHandler sh(myStrm);

  // attach it to logger
  int id = l.addHandler(sh);

  // get a stream and log something to it
  Lucee::LogStream dbgStrm = l.getDebugStream();
  dbgStrm << "Hello World";

  LC_ASSERT("Testing log message",
    myStrm.str() == "Hello World");

  // remove stream handler from logger
  LC_ASSERT("Testing if stream handler can be removed",
    l.removeHandler(id) == true);

  // now try to log to it
  dbgStrm << "This should not be logged";
  // stream should not have changed
  LC_ASSERT("Testing log message",
    myStrm.str() == "Hello World");
}

int
main(void)
{

  LC_BEGIN_TESTS("lclogger");
  test_a();
  LC_END_TESTS;
}
