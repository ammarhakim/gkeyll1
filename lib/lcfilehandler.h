/**
 * @file	lcfilehandler.h
 *
 * @brief	Class for loggers.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_FILE_HANDLER_H
#define LC_FILE_HANDLER_H

// lib includes
#include <lclogrecordhandler.h>

// std includes
#include <ios>
#include <fstream>

namespace Lucee
{

/**
 * Class to log messages to a file. By default the file is truncated
 * while opening it.
 */
  class FileHandler : public Lucee::LogRecordHandler
  {
    public:
/** 
 * Create a new file handler. This handler writes to a specified
 * output file.  In general, if the file exisits it will be truncated
 * and its existing contents lost.
 *
 * @param fname Name of the file.
 * @param mode File permissions.
 */
      FileHandler(const std::string& fname, std::ios_base::openmode mode=std::ios_base::trunc);

    private:
/** File stream for logging messages */
      std::ofstream _file;
  };
}

#endif // LC_LOG_RECORD_HANDLER_H
