/**
 * @file	LcFieldFactory.h
 *
 * @brief	A factory to make fields that live on grids.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_FIELD_FACTORY_H
#define LC_FIELD_FACTORY_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDataStructIfc.h>
#include <LcGenericFactory.h>

// std includes
#include <string>
#include <vector>

namespace Lucee
{
  template <unsigned NDIM>
  class FieldFactory : public Lucee::GenericFactory<Lucee::DataStructIfc>
  {
    public:
/** Class id: this is used by the registration system */
      static const char *id;

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      void readInput(Lucee::LuaTable& tbl);

/**
 * Create a new field that lives on a grid.
 *
 * @return pointer to new field.
 */
      Lucee::DataStructIfc* create(const Lucee::SolverIfc& solver);

    private:
/** Name of grid on which field lives */
      std::string onGrid;
/** Number of components in field */
      unsigned numComponents;
/** Ghost cells along lower side */
      int lowerGhost[NDIM];
/** Ghost cells along upper side */
      int upperGhost[NDIM];
  };
}

#endif // LC_FIELD_FACTORY_H
