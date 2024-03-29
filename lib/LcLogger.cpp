/**
 * @file	LcLogger.cpp
 *
 * @brief	Class for loggers.
 */

// lucee includes
#include <LcLogger.h>
#include <LcLogStream.h>

namespace Lucee
{
  Logger::Logger(const std::string &name) 
    : _parent(0), _name(name), _level(NOTSET), _oldLevel(NOTSET) 
  {
    _maxHandlerId = -1;
    // construct level -> string mapping
    _levelMap.insert( LevelPair_t(NOTSET, "notset") );
    _levelMap.insert( LevelPair_t(DEBUG, "debug") );
    _levelMap.insert( LevelPair_t(INFO, "info") );
    _levelMap.insert( LevelPair_t(WARNING, "warning") );
    _levelMap.insert( LevelPair_t(ERROR, "error") );
    _levelMap.insert( LevelPair_t(CRITICAL, "critical") );
    _levelMap.insert( LevelPair_t(DISABLED, "disabled") );

    // construct string -> level mapping
    _stringMap.insert( StringPair_t("notset", NOTSET) );
    _stringMap.insert( StringPair_t("debug", DEBUG) );
    _stringMap.insert( StringPair_t("info", INFO) );
    _stringMap.insert( StringPair_t("warning", WARNING) );
    _stringMap.insert( StringPair_t("error", ERROR) );
    _stringMap.insert( StringPair_t("critical", CRITICAL) );
    _stringMap.insert( StringPair_t("disabled", DISABLED) );
  }

  Logger::~Logger()
  {
    _handlers.erase(_handlers.begin(), _handlers.end());
  }

  std::string
  Logger::getName() const
  {
    return _name;
  }

  void 
  Logger::debug(const std::string& msg) const 
  {
    log(msg, DEBUG);
  }

  void 
  Logger::info(const std::string& msg) const 
  {
    log(msg, INFO);
  }

  void 
  Logger::warning(const std::string& msg) const 
  {
    log(msg, WARNING);
  }

  void 
  Logger::error(const std::string& msg) const 
  {
    log(msg, ERROR);
  }

  void
  Logger::critical(const std::string& msg) const 
  {
    log(msg, CRITICAL);
  }

  void 
  Logger::setLevel(LogMsgLevels level) 
  {
    _oldLevel = _level = level;
  }

  void 
  Logger::setLevel(const std::string& level) 
  {
    StringMap_t::const_iterator i =
      _stringMap.find(level);
    // set level, defaulting to NOTSET if incorrect string passed
    _oldLevel = _level = 
      ((i != _stringMap.end()) ? (*i).second : NOTSET);
  }

  LogMsgLevels 
  Logger::getLevel() const 
  {
    return _level;
  }

  std::string 
  Logger::getLevelStr() 
  {
    return _levelMap[_level];
  }

  void 
  Logger::disable() 
  {
    _oldLevel = _level;
    _level = DISABLED;
  }

  void 
  Logger::enable() 
  {
    _level = _oldLevel;
  }

  LogStream 
  Logger::getDebugStream() 
  {
    return LogStream(*this, DEBUG);
  }

  LogStream 
  Logger::getInfoStream()
  {
    return LogStream(*this, INFO);
  }

  LogStream 
  Logger::getWarningStream()
  {
    return LogStream(*this, WARNING);
  }

  LogStream 
  Logger::getErrorStream()
  {
    return LogStream(*this, ERROR);
  }

  LogStream 
  Logger::getCriticalStream()
  {
    return LogStream(*this, CRITICAL);
  }

  void
  Logger::log(const std::string& msg, LogMsgLevels withLevel) const 
  {
    // check if level of logger is sufficient to log this message
    if (_level <= withLevel)
    {
      // send message to each handler registered with logger
      std::map<int, Lucee::LogRecordHandler*>::const_iterator i;
      for (i=_handlers.begin(); i!=_handlers.end(); ++i) 
        i->second->write(msg);
    }
    // now send this same message to the parent if there is one
    if (_parent)
      _parent->log(msg, withLevel);
  }

  int
  Logger::addHandler(Lucee::LogRecordHandler& handler)
  {
    // add this to our list of handlers
    _handlers[++_maxHandlerId] = &handler;
    return _maxHandlerId;
  }

  bool
  Logger::removeHandler(int handlerId)
  {
    return _handlers.erase(handlerId) == 1;
  }
}
