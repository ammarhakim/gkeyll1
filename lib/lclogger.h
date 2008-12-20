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
#include <lclogrecordhandler.h>

namespace Lucee
{
// forward declare LogStream class
  class LogStream;

/**
 * Verbosity levels for log messages.
 */
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
/** Delcare LoggeBase as a friend so it access our privates */
      friend class Lucee::LoggerBase<Logger>;

/** Type for map of levels to strings */
      typedef std::map<unsigned, std::string, std::less<unsigned> > LevelMap_t;
/** Type for pairs of levels and strings */
      typedef std::pair<unsigned, std::string> LevelPair_t;

/** Type for map of strings to levels */
      typedef std::map<std::string, LogMsgLevels, std::less<std::string> > StringMap_t;
/** Type for pairs of strings and levels */
      typedef std::pair<std::string, LogMsgLevels>  StringPair_t;

/**
 * Creates a new logger with given name. The verbosity is set
 * NOTSET. However, by default no handlers are installed so no
 * messages will be logged. The name of the logger can be used as a
 * global handle to access this logger.
 *
 * @param name Name of logger.
 */
      Logger(const std::string &name);

/**
 * Log a debug message.
 *
 * @param msg Log message.
 */
      void debug(const std::string& msg) const;

/**
 * Log a info message
 *
 * @param msg Log message
 */
      void info(const std::string& msg) const;

/**
 * Log a warning message.
 *
 * @param msg Log message.
 */
      void warning(const std::string& msg) const;

/**
 * Log a error message.
 *
 * @param msg Log message.
 */
      void error(const std::string& msg) const;

/**
 * Log a critical message.
 *
 * @param msg Log message.
 */
      void critical(const std::string& msg) const;

/**
 * Generic messsage logger.
 *
 * @param msg Log message.
 * @param withLevel Level for this log message
 */
      void log(const std::string& msg, LogMsgLevels withLevel) const;

/**
 * Set verbosity level, passing an integer from Lucee::LogMsgLevels.
 *
 * @param level Level for this logger
 */
      void setLevel(LogMsgLevels level);

/**
 * Set verbosity level passing a string. Valid strings are: "notset",
 * "debug", "info", "warning", "error", "critical", "disabled".
 *
 * @param level Level string.
 */
      void setLevel(const std::string& level);

/**
 * Get verbosity level, one of Lucee::LogMsgLevels.
 *
 * @return verbosity level.
 */
      LogMsgLevels getLevel() const;

/**
 * Get verbosity level as a string.
 *
 * @return verbosity level.
 */
      std::string getLevelStr();

/**
 * Add a new handler to the logger.
 *
 * @param handler Instance of class LogRecordHandler. Note that the
 * handler is owned by the logger and hence is deleted when the logger
 * is deleted.
 */
      void addHandler(Lucee::LogRecordHandler& handler);

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
 * Return stream to log debug messages.
 *
 * @return log-stream for logging messages.
 */
      Lucee::LogStream getDebugStream();

/**
 * Return stream to log info messages.
 *
 * @return log-stream for logging messages.
 */
      Lucee::LogStream getInfoStream();

/**
 * Return stream to log warning messages.
 *
 * @return log-stream for logging messages.
 */
      Lucee::LogStream getWarningStream();

/**
 * Return stream to log error messages.
 *
 * @return log-stream for logging messages.
 */
      Lucee::LogStream getErrorStream();

/**
 * Return stream to log critical messages.
 *
 * @return log-stream for logging messages.
 */
      Lucee::LogStream getCriticalStream();

    private:
/**
 * Coping is not allowed.
 *
 * @param lg Logger to copy.
 */
      Logger(const Logger& lg);

/**
 * Assignment is not allowed.
 *
 * @param lg Logger to assign from.
 */
      Logger& operator=(const Logger& lg);

/** Pointer to parent logger */
      Logger *_parent;
/** Name of logger */
      std::string _name;
/** Level flags for logger verbosity */
      LogMsgLevels _level;
/** Previous level flags for logger verbosity */
      LogMsgLevels _oldLevel;
/** Map of levels to strings */
      LevelMap_t _levelMap;
/** Map of strings to level */
      StringMap_t _stringMap;
/** List of handlers regeistered with this logger */
      std::vector<Lucee::LogRecordHandler*> _handlers;
  };
}

#endif // LC_LOGGER_H
