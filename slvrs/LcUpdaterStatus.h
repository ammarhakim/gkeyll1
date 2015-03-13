/**
 * @file	LcUpdaterStatus.h
 *
 * @brief	Class indicating updater status.
 */

#ifndef LC_UPDATER_STATUS_H
#define LC_UPDATER_STATUS_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// gkeyll includes
#include <LcBasicObj.h>

// std includes
#include <string>

namespace Lucee
{
/**
 * Class to represent status of an updater's update() method. This
 * class has two methods, "getStatus" that returns true if the
 * update() method passes, or false otherwise. The second method is
 * getSuggestedDt() returns the time-step suggestion for the next
 * time-step.
 */
  class UpdaterStatus : public BasicObj
  {
    public:
/**
 * Create a new default updater status object: status is "true" and
 * suggestedDt is maximum allowable time-step.
 */
      UpdaterStatus();

/**
 * Create a new updater status object.
 *
 * @param status Status of update() method.
 * @param suggestedDt Suggested time-step.
 */
      UpdaterStatus(bool status, double suggestedDt);

/**
 * Create a new updater status object.
 *
 * @param status Status of update() method.
 * @param suggestedDt Suggested time-step.
 * @param msg Message for status.
 */
      UpdaterStatus(bool status, double suggestedDt,
        const std::string& msg);

/**
 * @return Status of update() method.
 */
      bool getStatus() const { return status; }

/**
 * @return suggested time-step for next time-step.
 */
      bool getSuggestedDt() const { return dt; }

/** Status of update */
      bool status;
/** Suggested time-step */
      double dt;
/** Message to explicate the status */
      std::string message;
  };
}

#endif // LC_UPDATER_STATUS_H
