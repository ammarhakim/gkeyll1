/**
 * @file	LcFileHandler.cpp
 *
 * @brief	Class for loggers.
 */

// lucee includes
#include <LcFileHandler.h>

// std includes
#include <fstream>

namespace Lucee
{
  FileHandler::FileHandler(const std::string& fname, std::ios_base::openmode mode)
    : Lucee::LogRecordHandler(_file), _file(fname.c_str(), mode)
  {
  }
}
