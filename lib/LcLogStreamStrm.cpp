/**
 * @file	LcLogStreamStrm.cpp
 *
 * @brief	Class for log streams.
 *
 * @version	$Id: LcLogStreamStrm.cpp 34 2008-12-20 23:58:04Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// lucee includes
#include <LcLogger.h>
#include <LcLogStream.h>

namespace Lucee
{
  LogStreamStrm::LogStreamStrm(Lucee::Logger& log, int level)
    : _logger(log), _level(level) 
  {
  }

  void
  LogStreamStrm::_logIt(const std::ostringstream& str)
  {
    _logger.log(str.str(), (Lucee::LogMsgLevels)_level);
  }

}
