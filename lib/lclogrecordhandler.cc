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
    // detach from all loggers
    std::vector<std::string> nms = loggerNames();
    for (unsigned i=0; i<nms.size(); ++i)
      detachFromLogger(nms[i]);
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
  LogRecordHandler::attachToLogger(const std::string& name) 
  {
    Lucee::Logger& logger = Lucee::Logger::get(name);
    attachToLogger(logger);
  }

  void
  LogRecordHandler::attachToLogger(Lucee::Logger& logger) 
  {
    int id = logger.addHandler(*this);
    _loggers[logger.getName()] = &logger;
    _handlerIds[logger.getName()] = id;
  }

  bool
  LogRecordHandler::detachFromLogger(const std::string& name) 
  {
    Lucee::Logger& logger = Lucee::Logger::get(name);
    return detachFromLogger(logger);
  }

  bool
  LogRecordHandler::detachFromLogger(Lucee::Logger& logger) 
  {
    std::string name = logger.getName();
    // check if this handler is really attached to the logger
    std::map<std::string, Logger*>::iterator itr
      = _loggers.find(name);
    if (itr != _loggers.end())
    {
      // detach handler from the logger
      itr->second->removeHandler(_handlerIds[name]);
      // remove entries from the map
      _loggers.erase(name);
      _handlerIds.erase(name);
      return true;
    }
    return false;
  }

  std::vector<std::string>
  LogRecordHandler::loggerNames() const
  {
    std::vector<std::string> nms;
    std::map<std::string, Logger*>::const_iterator itr;
    for (itr = _loggers.begin(); itr != _loggers.end(); ++itr)
      nms.push_back(itr->first);
    return nms;
  }
}
