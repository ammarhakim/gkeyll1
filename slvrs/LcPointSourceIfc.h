/**
 * @file	LcPointSourceIfc.h
 *
 * @brief	Interface to "point" sources.
 */
#ifndef LC_POINT_SOURCE_IFC_H
#define LC_POINT_SOURCE_IFC_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBasicObj.h>

// std includes
#include <vector>

namespace Lucee
{
/**
 * Base class for sources that depend on values in a single cell.
 */
  class PointSourceIfc : public Lucee::BasicObj
  {
    public:
/** 
 * Create new point source object that takes specified number of
 * inputs and computes specified number of outputs. If the flag
 * "allowArb" is set the number of input and output specified in the
 * constructor are discarded and set from the size of the
 * 'inpComponents' and 'outComponents' lists.
 */
      PointSourceIfc(unsigned nInp, unsigned nOut, bool allowArb = false);

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Get number of input variables.
 *
 * @param number of input variables.
 */
      unsigned getNumInput() const { return nInp; }

/**
 * Get number of output variables.
 *
 * @param number of output variables.
 */
      unsigned getNumOutput() const { return nOut; }

/**
 * Compute sources and store them in supplied output vector.
 *
 * @param loc Coordinate at which source is requested.
 * @param inp Input values at which source is requested.
 * @param src On output, source.
 */
      void calcSource(const double loc[3], const double *inp, double *src);

    protected:
/**
 * Compute sources and store them in supplied output vector. The
 * vector 'src' is pre-allocated. Derived class method should use the
 * getData() method to get data it needs in computing the sources.
 *
 * @param loc Coordinate at which source is requested.
 * @param src On output, source.
 */
      virtual void getSource(const double loc[3], std::vector<double>& src) = 0;

/**
 * Get data at i-th input component. This function takes care of
 * mapping between the global component list and the local component
 * list so that derived classes do not have to keep track of mapping
 * themselves.
 */
      double getData(unsigned idx) const
      {
        return data[inpComponents[idx]];
      }

    private:
/** Number of input variables */
      unsigned nInp;
/** Number of output variables */
      unsigned nOut;
/** Are arbitrary number of input/output allowed? */
      bool allowArb;
/** Pointer to data */
      const double *data;
/** List of input components */
      std::vector<unsigned> inpComponents;
/** List of output components */
      std::vector<unsigned> outComponents;
  };
}

#endif // LC_POINT_SOURCE_IFC_H
