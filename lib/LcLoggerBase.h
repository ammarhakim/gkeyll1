/**
 * @file	LcLoggerBase.h
 *
 * @brief	Base class for loggers.
 *
 * @version	$Id: LcLoggerBase.h 69 2008-12-30 07:35:17Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_LOGGER_BASE_H
#define LC_LOGGER_BASE_H

// lucee includes
#include <LcExcept.h>

// std includes
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace Lucee
{
/**
 * Base class for loggers. This class provides a method to register
 * and access loggers.
 */
  template<class T>
  class LoggerBase
  {
    public:
/**
 * Create a logger with a given name. A new logger is created and
 * returned. Note that calling this method with the same name will
 * always return the same logger object. Thus once loggers are created
 * they become global objects. The logger names are case sensitive.
 *
 * @param name name of logger to create. 
 * @return a new logger object.
 */    
      static T& create(const std::string& name);

/**
 * Returns a logger with a given name. If such a logger does not exist
 * an exception is thrown.
 *
 * @param name specifies the logger to get. 
 * @return the Logger object corresponding to the name.
 */
      static T& get(const std::string& name);

/**
 * Delete all loggers registered in the system.
 */
      static void cleanup();

    private:
/** Type for map of string to loggers */
      typedef std::map<std::string, T*, std::less<std::string> > LoggerMap_t;
/** Type for pair of string and loggers */
      typedef std::pair<std::string, T*> LoggerPair_t;

/** Map of logger names to loggers */
      static LoggerMap_t *loggers;
  };

  // initialize maps
  template<class T>
  std::map<std::string, T*, std::less<std::string> >* Lucee::LoggerBase<T>::loggers = NULL;

  template<class T>
  T&
  Lucee::LoggerBase<T>::create(const std::string& nm) 
  {
    if (!loggers)
    {
      loggers = new LoggerMap_t();
    }          

    // check if logger with given name is in map
    typename LoggerMap_t::iterator i = loggers->find(nm);
    if (i != loggers->end())
    {
      Lucee::Except ex("LoggerBase::create: Logger with name ");
      ex << nm << " already exists";
      throw ex;
    }

    // it is not, so start creation process

    // find location of '.'
    size_t len = 0;
    for (size_t i = nm.size(); i >= 0; --i, ++len)
    {
      if ( (nm[i] == '.') || (i == 0) )
      {
        // create logger with given name
        T *l = new T(nm);
        loggers->insert(LoggerPair_t(nm, l));
        // set its parent 
        if (i != 0)
        {
          l->_parent = &Lucee::LoggerBase<T>::get(nm.substr(0, nm.size()-len));
        }
        return *l;
      }
    }
  }

  template<class T>
  T&
  Lucee::LoggerBase<T>::get(const std::string& nm) 
  {
    if (!loggers)
    {
      loggers = new LoggerMap_t();
    }

    typename LoggerMap_t::iterator i = loggers->find(nm);
    if (i != loggers->end())
    {
      return *i->second;
    }
    // thow exception as logger not found
    Lucee::Except ex("LoggerBase::get: Logger with name ");
    ex << nm << " not found";
    throw ex;
  }

  template<class T>
  void 
  Lucee::LoggerBase<T>::cleanup()
  {
    typename LoggerMap_t::iterator itr;
    for (itr = loggers->begin(); itr != loggers->end(); ++itr)
    {
      delete itr->second;
    }
  }
}

#endif // LC_LOGGER_BASE_H
