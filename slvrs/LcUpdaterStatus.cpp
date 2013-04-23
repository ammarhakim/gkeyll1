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
    TxCommBase *comm = Loki::SingletonHolder<Lucee::Globals>
      ::Instance().comm;
  }

  UpdaterStatus::UpdaterStatus(bool myStatus, double mySuggestedDt)
  {
    TxCommBase *comm = Loki::SingletonHolder<Lucee::Globals>
      ::Instance().comm;
    int myLocalStatus = myStatus, myGlobalStatus;
    comm->allreduce(1, &myLocalStatus, &myGlobalStatus, TX_AND);
    status = myGlobalStatus;
    comm->allreduce(1, &mySuggestedDt, &dt, TX_MIN);
    //std::cout << comm->getRank() << " dt = " << mySuggestedDt << std::endl;

    if (status)
      message = "Updater status: Success";
    else
      message = "Updater status: Failure";
  }

  UpdaterStatus::UpdaterStatus(bool myStatus, double mySuggestedDt,
    const std::string& msg)
  {
    message = msg;

    TxCommBase *comm = Loki::SingletonHolder<Lucee::Globals>
      ::Instance().comm;
    int myLocalStatus = myStatus, myGlobalStatus;
    comm->allreduce(1, &myLocalStatus, &myGlobalStatus, TX_AND);
    status = myGlobalStatus;
    comm->allreduce(1, &mySuggestedDt, &dt, TX_MIN);
  }
}
