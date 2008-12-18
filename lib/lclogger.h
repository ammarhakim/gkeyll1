/**
 * @file	lclogger.h
 *
 * @brief	Class for loggers.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_LOGGER_H
#define LC_LOGGER_H

// lib includes
#include <lcexcept.h>
#include <lcloggerbase.h>

namespace Lucee
{
// forward declare LogStream class
  class Lucee::LogStream;

  enum LogMsgLevels 
  {
    NOTSET = 0,
    DEBUG = 10,
    INFO = 20,
    WARNING = 30,
    ERROR = 40,
    CRITICAL = 50,
    DISABLED = 1000
  };

  class Logger : public Lucee::LoggerBase<Logger>
  {
    public:
      friend class Lucee::LoggerBase<Logger>;

      // map type to map logger level to its string representation
      typedef std::map<unsigned, std::string, std::less<unsigned> > LevelMap_t;
      typedef std::pair<unsigned, std::string> LevelPair_t;

      typedef std::map<std::string, LogMsgLevels, std::less<std::string> > StringMap_t;
      typedef std::pair<std::string, LogMsgLevels>  StringPair_t;

/**
 * Creates a new logger with given name. The verbosity is set
 * NOTSET. However, by default no handlers are installed so no
 * messages will be logged.
 */
      Logger(const std::string &name);

/**
 * Delete logger and all log record handlers
 */
      virtual ~Logger();

/**
 * Log a debug message
 */
      void debug(const std::string& msg) const;

/**
 * Log a info message
 */
      void info(const std::string& msg) const;

/**
 * Log a warning message
 */
      void warning(const std::string& msg) const;

/**
 * Log a error message
 */
      void error(const std::string& msg) const;

/**
 * Log a critical message
 */
      void critical(const std::string& msg) const;

/**
 * Generic messsage logger.
 */
      void log(const std::string& msg, LogMsgLevels withLevel) const;

/**
 * Set verbosity level
 */
      void setLevel(LogMsgLevels level);

/**
 * Set verbosity level passing a string
 */
      void setLevel(const std::string& level);

/**
 * Get verbosity level
 */
      LogMsgLevels getLevel() const;

/**
 * Get verbosity level as a string
 */
      std::string getLevelStr();

/**
 * Add a new handler to the logger.
 *
 * @param handler Instance of class LogRecordHandler. Note that the
 * handler is not owned by the logger and hence is never deleted.
 */
      void addHandler(LogRecordHandler *handler);

/**
 * Disables all logging to this logger. Logging at the original level
 * can be resumed by calling the 'enable' method.
 */
      void disable();

/**
 * Re-enables logging to this logger. The verbosity level is set to
 * the one prior to the disable call.
 */
      void enable();

/**
 * Return stream to log debug messages
 */
      LogStream getDebugStream();

/**
 * Return stream to log info messages
 */
      LogStream getInfoStream();

/**
 * Return stream to log warning messages
 */
      LogStream getWarningStream();

/**
 * Return stream to log error messages
 */
      LogStream getErrorStream();

/**
 * Return stream to log critical messages
 */
      LogStream getCriticalStream();

    private:
      Logger *_parent; // parent logger
      std::string _name;
      LogMsgLevels _level, _oldLevel;
      LevelMap_t _levelMap;
      StringMap_t _stringMap;
      std::vector<LogRecordHandler*> _handlers;
  };
}

#endif // LC_LOGGER_H
