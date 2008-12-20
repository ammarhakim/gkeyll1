/**
 * @file	lcloggerbase.h
 *
 * @brief	Base class for loggers.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_LOGGER_BASE_H
#define LC_LOGGER_BASE_H

// lib includes
#include <lcexcept.h>

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
 * Returns a logger with a given name. 
 * 
 * If such a logger does not exist a new logger with that name is
 * created and returned. Note that calling this method with the same
 * name will always return the same logger object. Thus once loggers
 * are created they become global objects. The logger names are case
 * sensitive.
 *
 * @param name specifies the logger to get. 
 * @return the Logger object corresponding to the name.
 */    
      static T& get(const std::string& name);

/**
 * Returns the default logger.
 *
 * @return default logger
 */
      static T& getDefault();

/**
 * Returns a logger with a given name. If such a logger does not exist
 * an exception is thrown.
 *
 * @param name specifies the logger to get. 
 * @return the Logger object corresponding to the name.
 */
      static T& getSafe(const std::string& name);

/**
 * Delete all loggers registered in the system.
 */
      static void cleanUp();

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
  Lucee::LoggerBase<T>::get(const std::string& nm) 
  {
    // TODO, handle this using a const at class level.
    static const std::string defaultLoggerName("84BEF5E7-16AB-4FFF-A33B-36DD7140AA02");

    if (!loggers)
    {
      loggers = new LoggerMap_t();
      // Create default logger.
      T *l = new T(defaultLoggerName);
      loggers->insert(LoggerPair_t(defaultLoggerName, l));
    }          

    // check if logger with given name is in map
    typename LoggerMap_t::iterator i = loggers->find(nm);
    if (i != loggers->end())
    {
      return *i->second;
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
  Lucee::LoggerBase<T>::getDefault()
  {
    // TODO, handle this using a const at class level.
    static const std::string defaultLoggerName = "84BEF5E7-16AB-4FFF-A33B-36DD7140AA02";
    return *Lucee::LoggerBase<T>::get(defaultLoggerName);
  }

  template<class T>
  T&
  Lucee::LoggerBase<T>::getSafe(const std::string& nm) 
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
    Lucee::Except ex;
    ex << "Logger " << std::string(nm) << " not found";
    throw ex;
  }

  template<class T>
  void 
  Lucee::LoggerBase<T>::cleanUp()
  {
    typename LoggerMap_t::iterator itr;
    for (itr = loggers->begin(); itr != loggers->end(); ++itr)
    {
      delete itr->second;
    }
  }
}

#endif // LC_LOGGER_BASE_H
