/**
 * @file	LcOverlappingFieldCopy.h
 *
 * @brief	Updater to copy data into ghost cells between two overlapping fields
 */

#ifndef LC_OVERLAPPING_FIELD_COPY_H
#define LC_OVERLAPPING_FIELD_COPY_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcNodalFiniteElementIfc.h>
#include <LcPointSourceIfc.h>
#include <LcUpdaterIfc.h>

namespace Lucee
{
/**
 * Updater to copy data into ghost cells between two overlapping fields
 */
  template <unsigned NDIM>
  class OverlappingFieldCopy : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create new integrator.
 */
      OverlappingFieldCopy();

/** Dtor */
      ~OverlappingFieldCopy();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Advance the solution to specified time. Updaters that do not have a
 * concept of time should ignore the time parameter.
 *
 * @param t Time to advance the solution to.
 * @return Status of updater.
 */
      Lucee::UpdaterStatus update(double t);

/**
 * Declare the types of input and output variables accepted by this
 * updater. This must be provided by the derived classes for the
 * type-checking to pass. Inside the implementation of this method the
 * derived class must make a sequence of appendInpVarType() and
 * appendOutVarType() calls to declare the input/output data structure
 * types.
 */
      void declareTypes();

    private:
/* Number of cells that overlap */
      unsigned numCells;
/* Direction of overlap */
      unsigned dir;
  };
}

#endif // LC_OVERLAPPING_FIELD_COPY_H
