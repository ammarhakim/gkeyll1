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
  LC_ASSERT("Testing level", l.getLevelStr() == "debug");

  Lucee::LogStream ls = l.getDebugStream();
  ls << "Hello World" << std::endl;
}

int
main(void)
{
  // create a logger for testing
  Lucee::Logger& w = Lucee::Logger::get("lucee");
  w.setLevel(Lucee::DEBUG);
    
  Lucee::FileHandler fh("lucee-detailed.log");
  w.addHandler(fh);

  LC_BEGIN_TESTS("lclogger");
  test_a();
  LC_END_TESTS;
}
