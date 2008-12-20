/**
 * @file	lcstreamhandler.h
 *
 * @brief	Class for loggers.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_STREAM_HANDLER_H
#define LC_STREAM_HANDLER_H

// lib includes
#include <lclogrecordhandler.h>

// std includes
#include <ios>
#include <iostream>

namespace Lucee
{

/**
 * Class to log messages to an open stream.
 */
  class StreamHandler : public LogRecordHandler
  {
    public:
/**
 * Create new handler with the specified output stream.
 *
 * @param stream Stream to log messages to.
 */
      StreamHandler(std::ostream& stream);
  };

}

#endif // LC_STREAM_HANDLER_H
