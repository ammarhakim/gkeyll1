/**
 * @file	lclogstream.cc
 *
 * @brief	Class for log streams.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

// lib includes
#include <lclogger.h>
#include <lclogstream.h>

namespace Lucee
{
  LogStream::LogStream(Logger* log, int level)
    : useCount(new int(1))
  {
    strm = new LogStreamStrm(log, level);
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
