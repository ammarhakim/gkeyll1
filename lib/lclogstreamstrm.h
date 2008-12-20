/**
 * @file	lclogstreamstrm.h
 *
 * @brief	Class for log streams.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_LOG_STREAM_STRM_H
#define LC_LOG_STREAM_STRM_H

namespace Lucee
{
// forward declare Logger class
  class Logger;
// forward declare LogStream classes
  class LoggerStream;

  class LogStreamStrm 
  {
    public:
/** Let LogStream become our friend */
      friend class LogStream;

/**
 * Output supplied value
 *
 * @param val value to output
 * @return reference to this stream object
 */
      template <typename T>
      LogStreamStrm& operator<<(const T& val) 
      {
        _sstrm.str("");
        _sstrm << val;
        this->_logIt(_sstrm);
        return *this;
      }

/**
 * I/O for manipulators
 *
 * @param p manipulator object
 * @return reference to this stream object
 */
      LogStreamStrm&
      operator<<(std::ostream& (*p)(std::ostream&)) 
      {
        _sstrm.str("");
        _sstrm << p;
        this->_logIt(_sstrm);
        return *this;
      }

/**
 * I/O for manipulators
 *
 * @param p manipulator object
 * @return reference to this stream object
 */
      LogStreamStrm& 
      operator<<(std::ios& (*p)(std::ios&)) 
      {
        _sstrm.str("");
        _sstrm << p;
        this->_logIt(_sstrm);
        return *this;
      }

    private:
/** 
 * Ctor is private so only log-stream can make instances.
 */
      LogStreamStrm(Logger& log, int level);

/**
 * Copy ctor is private to avoid copying
 */
      LogStreamStrm(const LogStreamStrm&);

/** Output stream where messages are stored before being forwared to logstreams */
      std::ostringstream _sstrm;
/** Reference to logger*/
      Logger& _logger;
/** Level for log stream */
      int _level;

/** 
 * Log a message to the logger.
 *
 * @param str stream of log message.
 */
      void _logIt(const std::ostringstream& str);
  };
}

#endif // LC_LOG_STREAM_STRM_H
