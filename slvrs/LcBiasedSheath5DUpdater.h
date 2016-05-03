/**
 * @file	LcBiasedSheath5DUpdater.h
 *
 * @brief	Applies electrostatic logical sheath BCs to a 5D distribution function
 */

#ifndef LC_BIASED_SHEATH_5D_UPDATER
#define LC_BIASED_SHEATH_5D_UPDATER

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcMatrix.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

// eigen includes
#include <Eigen/Core>

namespace Lucee
{
/**
 * Applies particle refection BCs to distribution function
 */
  class BiasedSheath5DUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new projection updater */
      BiasedSheath5DUpdater();

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
/** Pointer to 5d nodal basis functions */
      Lucee::NodalFiniteElementIfc<5> *nodalBasis5d;
/** Pointer to 3d nodal basis functions */
      Lucee::NodalFiniteElementIfc<3> *nodalBasis3d;
/** Pointer to 2d nodal basis functions */
      Lucee::NodalFiniteElementIfc<2> *nodalBasis2d;
/** Mapping for 180 degree rotations */
      std::vector<unsigned> rotMap;
/** Contains the lower edge node numbers */
      std::vector<int> lowerEdgeNodeNums;
/** Contains the upper edge node numbers */
      std::vector<int> upperEdgeNodeNums;
/** Order of basis functions used */
      int polyOrder;
/** Keeps track of the offsets needed to get all nodes that share the same config. space location */
      std::vector<int> nodalStencil;
/** Used to compute first parallel velocity moment at a particular (x,y,z) */
      Eigen::MatrixXd momentMatrix;
/** Mass of species */
      double speciesMass;
/** Species charge */
      double speciesCharge;
/** Factor to multiply all results by (like 2*pi*B/m to account v_perp -> mu integration */
      double scaleFactor;
/** Stores location of gaussian quadrature points (on [-1,1]) for integration over entire 2d element */
      Eigen::MatrixXd gaussSurfCoords;
/** Stores gaussian quadrature weights for integration over 2d [-1,1]x[-1,1] element */
      Eigen::VectorXd gaussSurfWeights;

/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);

/**
 * Determines if two nodes have the same configuration space coordinates
 */
      bool sameConfigCoords(int srcIndex, int tarIndex, double dxMin, const Eigen::MatrixXd& nodeList);
  };
}

#endif // LC_BIASED_SHEATH_5D_UPDATER
