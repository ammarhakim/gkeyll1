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

WxLogStreamStrm::WxLogStreamStrm(WxLogger* log, int level)
  : _logger(log), _level(level) 
{
}

void
WxLogStreamStrm::_logIt(const std::wostringstream& str)
{
  _logger->log(str.str(), (WxLogger::eLevels)_level);
}
