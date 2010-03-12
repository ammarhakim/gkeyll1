/**
 * @file	LcLogStream.cpp
 *
 * @brief	Class for log streams.
 *
 * @version	$Id: LcLogStream.cpp 34 2008-12-20 23:58:04Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
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
