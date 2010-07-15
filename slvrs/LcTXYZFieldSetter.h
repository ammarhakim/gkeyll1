/**
 * @file	LcTXYZFieldSetter.h
 *
 * @brief	Linear combiner updater.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_TXYZ_FIELD_SETTER_H
#define LC_TXYZ_FIELD_SETTER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcFunctionIfc.h>
#include <LcUpdaterIfc.h>

namespace Lucee
{
/**
 * Updater to perform linear combination of input fields.
 */
  template <unsigned NDIM>
  class TXYZFieldSetter : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create new combiner.
 */
      TXYZFieldSetter();

/** Delete updater */
      ~TXYZFieldSetter();

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

/**
 * Set function object pointer to supplied one.
 *
 * @param f Function object pointer.
 */
      void setFunObj(Lucee::FunctionIfc& f);

    private:
/** Pointer to function object */
      Lucee::FunctionIfc *func;
/** Flag to indicate if this class owns func pointer */
      bool isOwner;
  };
}

#endif // LC_TXYZ_FIELD_SETTER_H
