/**
 * @file	LcLogStreamStrm.cpp
 *
 * @brief	Class for log streams.
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
