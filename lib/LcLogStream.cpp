/**
 * @file	LcLogStream.cpp
 *
 * @brief	Class for log streams.
 */

// lucee includes
#include <LcLogger.h>
#include <LcLogStream.h>

namespace Lucee
{
  LogStream::LogStream(Logger& log, int level)
    : useCount(new int(1)), strm(new LogStreamStrm(log, level))
  {
  }

  LogStream::LogStream(const LogStream& ls)
    : strm(ls.strm)
  {
    ++*ls.useCount;
    useCount = ls.useCount;
  }

  LogStream::~LogStream()
  {
    if (--*useCount == 0)
    {
      delete useCount;
      delete strm;
    }
  }

}
