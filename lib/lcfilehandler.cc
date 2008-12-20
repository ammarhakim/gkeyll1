/**
 * @file	lcfilehandler.cc
 *
 * @brief	Class for loggers.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

// lib includes
#include <lcfilehandler.h>
#include <fstream>

namespace Lucee
{
  FileHandler::FileHandler(const std::string& fname, std::ios_base::openmode mode)
    : Lucee::LogRecordHandler(_file), _file(fname.c_str(), mode)
  {
  }
}
