/**
 * @file	LcCopyBoundaryCondition.cpp
 *
 * @brief	Class for applying copy boundary conditions.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCopyBoundaryCondition.h>

namespace Lucee
{
// set costructor name
  const char *CopyBoundaryCondition::id = "Copy";

  void
  CopyBoundaryCondition::applyBc(
    const Lucee::RectCoordSys& c, const Lucee::ConstFieldPtr<double>& qin, Lucee::FieldPtr<double>& qbc)
  {
// just copy data over
    for (unsigned i=0; i<qin.getNumComponents(); ++i)
      qbc[i] = qin[i];
  }
}
