/**
 * @file	LcMusclHancock1DUpdater.h
 *
 * @brief	Solver for 1D Euler equations using MUSCl-Hancock scheme.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcMusclHancock1DUpdater.h>

namespace Lucee
{
/** Class id: this is used by registration system */
  const char *MusclHancock1DUpdater::id = "MusclHancock1D";

  MusclHancock1DUpdater::MusclHancock1DUpdater()
    : Lucee::UpdaterIfc() {
  }

  void
  MusclHancock1DUpdater::readInput(Lucee::LuaTable& tbl) {
// call base class method
    UpdaterIfc::readInput(tbl);
  }

  Lucee::UpdaterStatus
  MusclHancock1DUpdater::update(double t) {

    return Lucee::UpdaterStatus();
  }

  void
  MusclHancock1DUpdater::declareTypes() {
  }
}
