/**
 * @file	LcGlobals.h
 *
 * @brief	Class to hold global data used in various parts of Lucee.
 */

#ifndef LC_GLOBALS_H
#define LC_GLOBALS_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// txbase includes
#include <TxCommBase.h>

#ifdef HAVE_MPI
#include <TxMpiBase.h>
#else
#include <TxSelfBase.h>
#endif

// std includes
#include <string>

namespace Lucee
{
/**
 * Class to hold global data used in Lucee
 */
  struct Globals
  {
/** Create a new global object */
      Globals() 
      {
#ifdef HAVE_MPI
        comm = new TxMpiBase();
#else
        comm = new TxSelfBase();
#endif
      }

/** Delete global data */
      ~Globals()
      {
        delete comm;
      }

/** Output prefix for files */
      std::string outPrefix;
/** Pointer to communicator object */
      TxCommBase *comm;
  };
}

#endif // LC_GLOBALS_H
