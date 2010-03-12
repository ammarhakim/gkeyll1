/**
 * @file	LcFileHandler.cpp
 *
 * @brief	Class for loggers.
 *
 * @version	$Id: LcFileHandler.cpp 34 2008-12-20 23:58:04Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
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
