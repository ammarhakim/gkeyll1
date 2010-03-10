/**
 * @file	LcRteHomogeneousSlab.h
 *
 * @brief	Radiative transfer equation in homogeneous slab.
 *
 * @version	$Id: LcRteHomogeneousSlab.h 331 2010-03-10 03:55:02Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_RTE_HOMOGENEOUS_SLAB_H
#define LC_RTE_HOMOGENEOUS_SLAB_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSolverIfc.h>

// std includes
#include <string>

namespace Lucee
{
/**
 * Class to solve the radiative transfer equation in a homogeneous slab.
 */
  class RteHomogeneousSlab : public Lucee::SolverIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create a new RET solver object.
 */
      RteHomogeneousSlab();

/**
 * Destroy solver object.
 */
      virtual ~RteHomogeneousSlab();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Bootstrap method: Allocate data for solver.
 */
      virtual void buildData();

/**
 * Initialize algorithms needed for solver.
 */
      virtual void buildAlgorithms();

/**
 * Initialize solver, i.e. setup initial conditions. At the end of
 * this call, the solver should be ready for evolving the solution.
 */
      virtual void initialize();

/**
 * Advance the solution to specified time. Solvers that do not have a
 * concept of time should ignore the time parameter.
 *
 * @param t Time to advance the solution to.
 * @return Status of solution.
 */
      virtual int advance(double t);

/**
 * Write solver data to file.
 *
 * @param baseName Base name of output files. This should serve as a
 *   prefix for all output files.
 */
      virtual void writeToFile(const std::string& baseName) const;

/**
 * Restore solver data from file. This is called instead of the
 * initialize() method if the simulation is being restarted.
 *
 * @param baseName Base name of input files. This should serves as a
 *   prefix for all input files.
 */
      virtual void restoreFromFile(const std::string& baseName);

/**
 * Finalize solver: free resources, deallocate memory, close files
 * etc.
 */
      virtual void finalize();
    private:
/** Number of expansion coefficients */
      unsigned L;
/** Number of quadrature points in each hemisphere */
      unsigned N;
/** Cosine of incidence angle */
      double mu0;
/** Optical depth */
      double tau0;
  };
}

#endif // LC_RTE_HOMOGENEOUS_SLAB_H
