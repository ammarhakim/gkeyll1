/**
 * @file	LcStreamHandler.cpp
 *
 * @brief	Class for loggers.
 *
 * @version	$Id: LcStreamHandler.cpp 34 2008-12-20 23:58:04Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// lucee includes
#include <LcStreamHandler.h>

namespace Lucee
{
  StreamHandler::StreamHandler(std::ostream& stream)
    : Lucee::LogRecordHandler(stream)
  {
  }
}
