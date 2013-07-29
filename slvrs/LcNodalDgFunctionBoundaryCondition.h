/**
 * @file	LcNodalDgFunctionBoundaryCondition.h
 *
 * @brief	Class for applying boundary conditions specified via a Lua function.
 */

#ifndef LC_NODAL_DG_FUNCTION_BOUNDARY_CONDITION_H
#define LC_NODAL_DG_FUNCTION_BOUNDARY_CONDITION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBoundaryCondition.h>
#include <LcMatrix.h>
#include <LcNodalFiniteElementIfc.h>

namespace Lucee
{
/** Interface class for applying boundary conditions */
  template <unsigned NDIM>
  class NodalDgFunctionBoundaryCondition : public Lucee::BoundaryCondition
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/** 
 * Apply boundary conditions. The variable 'qin' represents data in
 * the skin cell (first/last interior cell) while 'qbc' represents
 * data in the ghost cells.
 *
 * @param tm Time at which source is requested.
 * @param loc Coordinate at which source is requested.
 * @param idx Grid index at which source is requested.
 * @param c Coordinate system to use.
 * @param qin Data in skin cell.
 * @param qbc Data in ghost cell.
 */
      void applyBc(double tm, const double loc[3], const int *idx,
        const Lucee::RectCoordSys& c, 
        const Lucee::ConstFieldPtr<double>& qin, Lucee::FieldPtr<double>& qbc);

    private:
/** Reference to Lua function */
      int fnRef;
/** Nodal finite element basis functions */
      Lucee::NodalFiniteElementIfc<NDIM> *nodalBasis;
/** Nodal coordinates */
      Lucee::Matrix<double> nodeCoords;
/** Unit mapping (I am not sure what this is for) */
      std::vector<int> ndIds;

/**
 * Evaluate function at specified location and fill output array with
 * result.
 *
 * @param L Lua state object to use.
 * @param tm Time to evaluate function at.
 * @param nc Matrix with nodal cooridates.
 * @param nn Node number
 * @param res On output, result of evaluating function.
 */
      void evaluateFunction(Lucee::LuaState& L, double tm, 
        const Lucee::Matrix<double>& nc, unsigned nn, std::vector<double>& res);
  };
}

#endif // LC_NODAL_DG_FUNCTION_BOUNDARY_CONDITION_H
