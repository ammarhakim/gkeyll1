/**
 * @file	lcstreamhandler.cc
 *
 * @brief	Class for loggers.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

// lib includes
#include <lcstreamhandler.h>

namespace Lucee
{
  StreamHandler::StreamHandler(std::ostream& stream)
    : Lucee::LogRecordHandler(stream)
  {
  }
}
