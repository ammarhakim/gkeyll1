/**
 * @file	lclogger.cxx
 *
 * @brief	Unit tests for logging classes
 */

// lucee include
#include <LcFileHandler.h>
#include <LcLogStream.h>
#include <LcLogger.h>
#include <LcStreamHandler.h>
#include <LcTest.h>

void
test_a()
{
  Lucee::Logger& l = Lucee::Logger::create("lucee");

// set level for logger
  l.setLevel("debug");
  LC_ASSERT("Testing level", l.getLevelStr() == "debug");

// create a stream handler
  std::ostringstream myStrm;
  Lucee::StreamHandler strmHandler(myStrm);
  strmHandler.attachToLogger("lucee");

// get a stream and log something to it
  Lucee::LogStream dbgStrm = l.getDebugStream();
  dbgStrm << "Hello World";

  LC_ASSERT("Testing log message",
    myStrm.str() == "Hello World");

// remove stream handler fromm logger
  LC_ASSERT("Checking if detaching from logger works",
    strmHandler.detachFromLogger("lucee") == true);

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

