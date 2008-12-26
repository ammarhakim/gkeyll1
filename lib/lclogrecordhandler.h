/**
 * @file	lclogrecordhandler.h
 *
 * @brief	Class for loggers.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_LOG_RECORD_HANDLER_H
#define LC_LOG_RECORD_HANDLER_H

// std includes
#include <ios>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace Lucee
{
// forward declare logger class
  class Logger;

/**
 * Class to handle log messages generated by the logging system. This
 * class can not be used directly. Instead, children handlers must
 * provide the needed functionality.
 */
  class LogRecordHandler
  {
    public:
/**
 * Destroy handler.
 */
      virtual ~LogRecordHandler();

/**
 * Write message 'msg' to the output stream.
 *
 * @param msg Message to write to stream.
 */
      virtual void write(const std::string& msg);

/**
 * Attach handler to a logger.
 *
 * @param name name of logger to which this handler should be attached.
 */
      void attachToLogger(const std::string& name);

/**
 * Attach handler to a logger.
 *
 * @param logger logger to which this handler should be attached.
 */
      void attachToLogger(Lucee::Logger& logger);

/**
 * Detach from logger.
 *
 * @param logger name of logger from which this handler should be detached.
 * @return true if detach worked, false otherwise.
 */
      bool detachFromLogger(const std::string& name);

/**
 * Detach from logger.
 *
 * @param logger logger from which this handler should be detached.
 * @return true if detach worked, false otherwise.
 */
      bool detachFromLogger(Lucee::Logger& logger);

/**
 * Get names of loggers this handler is attached to.
 *
 * @return names of loggers this handler is attached to.
 */
      std::vector<std::string> loggerNames() const;

    protected:
/**
 * Create new handler, attaching it to the specified stream.
 *
 * @param stream Stream to attach to handler.
 */
      LogRecordHandler(std::ostream& stream);

    private:
/** Reference to output stream */
      std::ostream& _outStream;
/** Map of logger names to loggers */
      std::map<std::string, Logger*> _loggers;
/** Map of logger names to handler ids */
      std::map<std::string, int> _handlerIds;
  };
}

#endif // LC_LOG_RECORD_HANDLER_H
