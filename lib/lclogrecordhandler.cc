/**
 * @file	lclogrecordhandler.cc
 *
 * @brief	Class for loggers.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

// lib includes
#include <lclogger.h>
#include <lclogrecordhandler.h>

namespace Lucee
{
  LogRecordHandler::~LogRecordHandler()
  {
  }

  LogRecordHandler::LogRecordHandler(std::ostream& stream)
    : _outStream(stream)
  {
  }

  void
  LogRecordHandler::write(const std::string& msg) 
  {
    _outStream << msg;
    _outStream << std::flush;
  }

  void
  LogRecordHandler::attachLogger(Lucee::Logger& logger) 
  {
    _loggers[logger.getName()] = &logger;
  }

  unsigned
  LogRecordHandler::numAttached() const
  {
    return _loggers.size();
  }
}
