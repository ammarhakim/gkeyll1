/**
 * @file	LcStreamHandler.cpp
 *
 * @brief	Class for loggers.
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
