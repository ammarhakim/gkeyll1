/**
 * @file	LcUpdaterStatus.cpp
 *
 * @brief	Class indicating updater status.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcUpdaterStatus.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <limits>

namespace Lucee
{
  UpdaterStatus::UpdaterStatus()
    : status(true), dt(std::numeric_limits<double>::max())
  {
    message = "Updater status: Success";
  }

  UpdaterStatus::UpdaterStatus(bool myStatus, double mySuggestedDt)
    : status(myStatus), dt(mySuggestedDt)
  {
  }

  UpdaterStatus::UpdaterStatus(bool myStatus, double mySuggestedDt,
    const std::string& msg)
    : status(myStatus), dt(mySuggestedDt)
  {
  }
}
