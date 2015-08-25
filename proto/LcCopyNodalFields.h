/**
 * @file	LcNodalVlasovUpdater.h
 *
 * @brief	Updater to copy one nodal field to another
 */

#ifndef LC_COPY_NODAL_FIELDS_H
#define LC_COPY_NODAL_FIELDS_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCDIM.h>
#include <LcField.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

namespace Lucee
{
/**
 * Copy one nodal field to another. The fields can be of different
 * dimensions. Input file is of dimension CDIM and output field of
 * dimension PDIM.
 */
  template <unsigned CDIM, unsigned PDIM>
  class CopyNodalFieldsUpdater : public Lucee::UpdaterIfc
  {
// Number of components for coordinate arrays etc.
      static const unsigned PNC = Lucee::CDIM<PDIM>::N;
      static const unsigned CNC = Lucee::CDIM<CDIM>::N;
      
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new nodal DG solver */
      CopyNodalFieldsUpdater();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Initialize solver, i.e. setup initial conditions. At the end of
 * this call, the solver should be ready for evolving the solution.
 */
      virtual void initialize();

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
/** Pointer to phase-space basis functions to use */
      Lucee::NodalFiniteElementIfc<PDIM> *phaseBasis;
/** Pointer to configuration-space basis functions to use */
      Lucee::NodalFiniteElementIfc<CDIM> *confBasis;
/** Mapping of node in phase-space to node in configuration space */
      std::vector<unsigned> phaseConfMap;      

/**
 * Check if coordinates of phase-space node and configuration-space
 * node are the same.
 */
      bool sameConfigCoords(unsigned n, unsigned cn, double dxMin,
        const Lucee::Matrix<double>& phaseC, const Lucee::Matrix<double>& confC);      
  };
}

#endif // LC_COPY_NODAL_FIELDS_H
