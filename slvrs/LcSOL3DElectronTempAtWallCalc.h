/**
 * @file  LcSOL3DElectronTempAtWallCalc.h
 *
 * @brief	Used to compute electron temperature at wall when sheath boundary
 * conditions are used by taking into account potential drop and cutoff velocity
 */

#ifndef LC_SOL_3D_ELECTRON_TEMP_AT_WALL_CALC
#define LC_SOL_3D_ELECTRON_TEMP_AT_WALL_CALC

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcMatrix.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>
#include <LcVector.h>
#include <LcDynVector.h>

// eigen includes
#include <Eigen/Dense>
#include <Eigen/LU>

namespace Lucee
{
/**
 * Applies particle refection BCs to distribution function
 */
  class SOL3DElectronTempAtWallCalc : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new projection updater */
      SOL3DElectronTempAtWallCalc();

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
/** Pointer to nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<3> *nodalBasis;
/** Pointer to 1D nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<1> *nodalBasis1d;
/** Electron mass */
      double elcMass;
/** Conversion of joules to eV */
      double eV;
/** Background magnetic field */
      double B0;
/** Contains the right edge node numbers */
      std::vector<int> rightEdgeNodeNums;
/** Contains the left edge node numbers */
      std::vector<int> leftEdgeNodeNums;
/**
 * Matrix of surface gaussian quadrature locations on bottom face..
 * There are three columns by default for (x,y,z)
 * and each row is a different quadrature point for doing surface integrals.
 */
      Eigen::MatrixXd gaussEdgeOrdinates;
/** Weights for edge gaussian quadrature points */
      std::vector<double> gaussEdgeWeights;
/** 
 * Interpolation matrix for bringing data that lives on the left or right edge
 * to quadrature points on same surface. It is constructed from the right edge
 * basis function evaluations but should work for the left edge as well.
 * The quadrature locations for these nodes are also in gaussSurfOrdinates and
 * gaussSurfWeights.
 */
      Eigen::MatrixXd edgeNodeInterpMatrix;

/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);
  };
}

#endif // LC_SOL_3D_ELECTRON_TEMP_AT_WALL_CALC
