/**
 * @file	LcStreamHandler.h
 *
 * @brief	Class for loggers.
 *
 * @version	$Id: LcStreamHandler.h 34 2008-12-20 23:58:04Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_STREAM_HANDLER_H
#define LC_STREAM_HANDLER_H

// lucee includes
#include <LcLogRecordHandler.h>

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
