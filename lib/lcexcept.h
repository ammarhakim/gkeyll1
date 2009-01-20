/**
 * @file	lcexcept.h
 *
 * @brief	Class to represent exceptions in Lucee
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_EXCEPT_H
#define LC_EXCEPT_H

// std includes
#include <sstream>
#include <exception>
#include <string>

namespace Lucee
{
/**
 * Lucee::Except is the class to use for creating and throwing
 * exceptions. It supports the standard I/O library operators to
 * create the exception messages. For example, one can do
 *
 * Lucee:Except ex("An exception ");
 * ex << "occured at line " << 25 << std::endl;
 */
  class Except : public std::exception
  {
    public:
/**
 * Create a new exception object with empty message.
 */
      Except();

/**
 * Create new exception with supplied message.
 *
 * @param str exception message.
 */
      Except(const std::string& str);

/**
 * Copy expection from supplied expection object.
 *
 * @param ex Excpetion to copy from
 */
      Except(const Except& ex);

/**
 * Destroy exception.
 */
      virtual ~Except() throw();

/**
 * Assignment operator from supplied expection object.
 *
 * @param ex Excpetion to copy from
 */
      Except& operator=(const Except& ex);

/**
 * Return error message: this is provided for compatibility with
 * std::exception interface.
 *
 * @return Message from exception.
 */
      virtual const char* what() const throw();

      template <typename T>
      Except& operator<<(const T& ex)
      {
        exceptStrm << ex;
        return *this;
      }

      Except& operator<<(std::ostream& (*p)(std::ostream&))
      {
        exceptStrm << p;
        return *this;
      }

    private:
/** Stream in which message is stored */
      mutable std::ostringstream exceptStrm;
/** String representing exception */
      mutable std::string exceptMsg;
  };
}

#endif // LC_EXCEPT_H
