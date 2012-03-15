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

// std includes
#include <limits>

namespace Lucee
{
  UpdaterStatus::UpdaterStatus()
    : status(true), dt(std::numeric_limits<double>::max())
  {
    message = "Updater status: Success";
  }

  UpdaterStatus::UpdaterStatus(bool status, double suggestedDt)
    : status(status), dt(suggestedDt)
  {
    if (status)
      message = "Updater status: Success";
    else
      message = "Updater status: Failure";
  }

  UpdaterStatus::UpdaterStatus(bool status, double suggestedDt,
    const std::string& msg)
    : status(status), dt(suggestedDt), message(msg)
  {
  }
}
