/**
 * @file	lclogstreamstrm.h
 *
 * @brief	Class for log streams.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_LOG_STREAM_STRM_H
#define LC_LOG_STREAM_STRM_H

namespace Lucee
{
// forward declare Logger class
  class Lucee::LoggerStream;
  class Lucee::Logger;

  class LogStreamStrm 
  {
    public:
      friend class Lucee::LogStream;

/**
 * Output supplied value
 *
 * @param val value to output
 * @return reference to this stream object
 */
      template <typename T>
      LogStreamStrm& operator<<(T val) 
      {
        _sstrm.str(L"");
        _sstrm << val;
        this->_logIt(_sstrm);
        return *this;
      }

/**
 * I/O for manipulators
 *
 * @param p manipulator object
 * @return reference to this stream object
 */
      LogStreamStrm&
      operator<<(std::ostream& (*p)(std::ostream&)) 
      {
        _sstrm.str(L"");
        _sstrm << p;
        this->_logIt(_sstrm);
        return *this;
      }

/**
 * I/O for manipulators
 *
 * @param p manipulator object
 * @return reference to this stream object
 */
      LogStreamStrm& 
      operator<<(std::ios& (*p)(std::ios&)) 
      {
        _sstrm.str(L"");
        _sstrm << p;
        this->_logIt(_sstrm);
        return *this;
      }

    private:
/** 
 * Ctor is private so only log-stream can make instances.
 */
      LogStreamStrm(Logger* log, int level);

/**
 * Copy ctor is private to avoid copying
 */
      LogStreamStrm(const LogStreamStrm&);

      std::wostringstream _sstrm;
      Logger *_logger;
      int _level;

      void _logIt(const std::wostringstream& str);
  };
}

#endif // LC_LOG_STREAM_STRM_H
