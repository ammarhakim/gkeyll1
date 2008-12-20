/**
 * @file	lclogstreamstrm.cc
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
