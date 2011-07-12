/**
 * @file	LcLogStream.h
 *
 * @brief	Class for loggers.
 */

#ifndef LC_LOG_STREAM_H
#define LC_LOG_STREAM_H

// lucee includes
#include <LcLogStreamStrm.h>

// std includes
#include <sstream>
#include <iostream>

namespace Lucee
{
// forward declare logger class
  class Lucee::Logger;

/**
 * Class that represents a output stream for logging.
 */
  class LogStream 
  {
    public:
/** 
 * Create new stream object.
 *
 * @param log reference to logger to which this stream is attached.
 * @param level Verbosity level of logger.
 */
      LogStream(Logger& log, int level);

/** 
 * Copy log stream from another logstream
 *
 * @param ls Logstream object to copy.
 */
      LogStream(const LogStream& ls);

/** 
 * Delete stream 
 */
      ~LogStream();

/**
 * Output supplied value
 *
 * @param val value to output
 * @return reference to this stream object
 */
      template <typename T>
      LogStream& operator<<(T val) 
      {
        strm->operator<<(val);
        return *this;
      }

/**
 * I/O for manipulators
 *
 * @param p manipulator object
 * @return reference to this stream object
 */
      LogStream& 
      operator<<(std::ostream& (*p)(std::ostream&))
      {
        strm->operator<<(p);
        return *this;
      }

/**
 * I/O for manipulators
 *
 * @param p manipulator object
 * @return reference to this stream object
 */
      LogStream& 
      operator<<(std::ios& (*p)(std::ios&)) 
      {
        strm->operator<<(p);
        return *this;
      }

    private:
/** For reference counting */
      mutable int *useCount;
/** Pointer to stream which does output */
      LogStreamStrm* strm;
  };
}

#endif // LC_LOG_STREAM
