/**
 * @file	LcSerendipityElement2D.h
 *
 * @brief       Reference finite element with serendipity basis
 */

#ifndef LC_SERENDEPITY_ELEMENT_2D_H
#define LC_SERENDEPITY_ELEMENT_2D_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcNodalFiniteElementIfc.h>
#include <LcRowMajorIndexer.h>

namespace Lucee
{
/**
 * Serendipity element in 2D. The reference element is a square [-1,1]
 * X [-1,1]. The supported elements are shown below with the
 * corresponding node layout.
 *
 * Order 1 element
 *
 *             4         3
 *             o----------o
 *             |          |
 *             |          |
 *             |          |
 *             |          |
 *             o----------o
 *             1          2
 *
 * Order 2 element
 *
 *             4     7     3
 *             o-----o-----o
 *             |           |
 *             |           |
 *           8 o           o 6
 *             |           |
 *             |           |
 *             o-----o-----o
 *             1     5     2
 */
  class SerendipityElement2D : public Lucee::NodalFiniteElementIfc<2>
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create a new serendipity element. This does not create a usable
 * object which can only be created from Lua.
 */
      SerendipityElement2D();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Get number of surface nodes along lower face in specified
 * direction.
 *
 * @param dir Direction to which face is perpendicular.
 * @return number of nodes in element.
 */
      unsigned getNumSurfLowerNodes(unsigned dir) const;

/**
 * Get number of surface nodes along upper face in specified
 * direction.
 *
 * @param dir Direction to which face is perpendicular.
 * @return number of nodes in element.
 */
      unsigned getNumSurfUpperNodes(unsigned dir) const;

/**
 * Get number of global nodes in element.
 *
 * @return number of nodes in element.
 */
      unsigned getNumGlobalNodes() const;

/**
 * Get local indices of nodes exclusively owned by each cell.
 *
 * @param ndIds On output indices. Vector is cleared and data filled in.
 */
      void getExclusiveNodeIndices(std::vector<unsigned>& ndIds);

/**
 * Get mapping of local node numbers in the current cell to global
 * node number. The input vector must be pre-allocated.
 *
 * @param lgMap Local node number to global node number mapping.
 */
      void getLocalToGlobal(std::vector<int>& lgMap) const;

/**
 * Get mapping of local node numbers to global node numbers on lower
 * face of element in direction 'dir' in the current cell . The input
 * vector must be pre-allocated.
 *
 * @param dir Direction to which face is perpendicular.
 * @param lgMap Local node number to global node number mapping.
 */
      void getSurfLowerLocalToGlobal(unsigned dir,
        std::vector<int>& lgMap) const;

/**
 * Get mapping of local node numbers to global node numbers on upper
 * face of element in direction 'dim' in the current cell . The input
 * vector must be pre-allocated.
 *
 * @param lgMap Local node number to global node number mapping.
 */
      void getSurfUpperLocalToGlobal(unsigned dim,
        std::vector<int>& lgMap) const;

/**
 * Get node numbers of the nodes on specified face of element. The
 * output vector must be pre-allocated.
 *
 * @param dir Direction to which face is perpendicular.
 * @param nodeNum Node numbers on face.
 */
      void getSurfLowerNodeNums(unsigned dir,
        std::vector<int>& nodeNum) const;

/**
 * Get node numbers of the nodes on specified face of element. The
 * output vector must be pre-allocated.
 *
 * @param dir Direction to which face is perpendicular.
 * @param nodeNum Node numbers on face.
 */
      void getSurfUpperNodeNums(unsigned dir,
        std::vector<int>& nodeNum) const;

/**
 * Get coordinates of all nodes in element. The output matrix
 * 'nodeCoords' should be pre-allocated have shape numNodes X 3.
 *
 * @param nodeCoords Node coordinates. Should be pre-allocated.
 */
      void getNodalCoordinates(Lucee::Matrix<double>& nodeCoords);

/**
 * Get weights for quadrature. The output vector should be
 * pre-allocated.
 *
 * @param w Weights for quadrature.
 */
      void getWeights(std::vector<double>& w);

/**
 * Get weights for quadrature on upper face. The output vector should
 * be pre-allocated.
 *
 * @param dir Direction to which face is perpendicular.
 * @param w Weights for quadrature.
 */
      void getSurfUpperWeights(unsigned dir,
        std::vector<double>& w);

/**
 * Get weights for quadrature on upper face. The output vector should
 * be pre-allocated.
 *
 * @param dir Direction to which face is perpendicular.
 * @param w Weights for quadrature.
 */
      void getSurfLowerWeights(unsigned dir,
        std::vector<double>& w);

/**
 * Get mass matrix for this reference element. The output matrix
 * should be pre-allocated.
 *
 * @param NjNk On output, mass matrix of element.
 */
      void getMassMatrix(Lucee::Matrix<double>& NjNk) const;

/**
 * Get mass matrix for this reference element. The output matrix
 * should be pre-allocated.
 *
 * @param dir Direction of face.
 * @param NjNk On output, mass matrix of element.
 */
      void getLowerFaceMassMatrix(unsigned dir,
        Lucee::Matrix<double>& NjNk) const;

/**
 * Get mass matrix for this reference element. The output matrix
 * should be pre-allocated.
 *
 * @param dir Direction of face.
 * @param NjNk On output, mass matrix of element.
 */
      void getUpperFaceMassMatrix(unsigned dir,
        Lucee::Matrix<double>& NjNk) const;

/**
 * Get stiffness matrix (grad.Nj \dot grad.Nk) for this reference
 * element. The output matrix should be pre-allocated.
 *
 * @param DNjDNk On output, stiffness matrix of element.
 */
      void getStiffnessMatrix(Lucee::Matrix<double>& DNjDNk) const;

/**
 * Get partial stiffness matrix (grad.Nj Nk) for this reference
 * element. The output matrix should be pre-allocated.
 *
 * @param dir Direction for gradient.
 * @param DNjNk On output, partial stiffness matrix of element.
 */
      void getGradStiffnessMatrix(unsigned dir, Lucee::Matrix<double>& DNjNk) const;

/**
 * Get data needed for Gaussian quadrature of specified order on this
 * element. All output matrices and vectors must be pre-allocated.
 *
 * @param nord Number of nodes in each direction.
 * @param interpMat On output, interpolation matrix.
 * @param ordinates On output, quadrature ordinates.
 * @param weights On output, quadrature weights.
 */
      void getGaussQuadData(unsigned norder, Lucee::Matrix<double>& interpMat,
        Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const;

/**
 * Extract nodal data at current grid location from field and copy it
 * into a vector. This basically "flattens" the nodal data consistent
 * with the node layout and the stiffness, mass matrices. The output
 * vector should be pre-allocated.
 *
 * @param fld Field to extract nodal data from.
 * @param data On output, this containts a copy of extracted data.
 */
      void extractFromField(const Lucee::Field<2, double>& fld,
        std::vector<double>& data);

/**
 * Copy all nodal data from field and put it into the data array. The
 * data pointer should be pre-allocated.
 *
 * @param fld Field to extract data from.
 * @param data Data space to copy into.
 */
      void copyAllDataFromField(const Lucee::Field<2, double>& fld, double *data);

/**
 * Copy all nodal data to field from a data array.
 *
 * @param data Data space to copy from.
 * @param fld Field to copy data to.
 */
      virtual void copyAllDataToField(const double *data, Lucee::Field<2, double>& fld);

    private:
/** Order of polynomial in element */
      unsigned polyOrder;
/** Mass matrix in reference coordinates */
      Lucee::Matrix<double> refNjNk;
/** Face-mass matrix in reference coordinates */
      Lucee::Matrix<double> refFaceNjNk_xl;
/** Face-mass matrix in reference coordinates */
      Lucee::Matrix<double> refFaceNjNk_xu;
/** Face-mass matrix in reference coordinates */
      Lucee::Matrix<double> refFaceNjNk_yl;
/** Face-mass matrix in reference coordinates */
      Lucee::Matrix<double> refFaceNjNk_yu;
/** Stiffness matrix in reference coordinates */
      Lucee::Matrix<double> refDNjDNk;
/** Grad-Stiffness (in direction X) matrix in reference coordinates */
      Lucee::Matrix<double> refDNjNk_0;
/** Grad-Stiffness (in direction Y) matrix in reference coordinates */
      Lucee::Matrix<double> refDNjNk_1;
/** Total global degrees of freedom */
      unsigned numGlobal;
/** Indexer to map (i,j) into a linear index */
      Lucee::RowMajorIndexer<2> idxr;
/** Number of cells in X and Y direction */
      unsigned numX, numY;
/** Weights for quadrature */
      std::vector<double> weights;
/** Weights for surface quadrature on 0-direction surfaces */
      std::vector<double> surfWeightsDir0;
/** Weights for surface quadrature on 1-direction surfaces */
      std::vector<double> surfWeightsDir1;

/**
 * Struct to hold data for Guassian quadrature.
 */
      struct GaussQuadData
      {
/**
 * Create object.
 * 
 * @param nord Numer of ordinates in each direction.
 * @param nlocal Total number of local nodes.
 */
          GaussQuadData(unsigned nord, unsigned nlocal)
            : ords(nord*nord, 3), weights(nord*nord), interpMat(nord*nord, nlocal)
          {
          }

/**
 * Reset object.
 * 
 * @param nord Numer of ordinates in each direction.
 * @param nlocal Total number of local nodes.
 */
          void reset(unsigned nord, unsigned nlocal)
          {
            ords = Lucee::Matrix<double>(nord*nord, 3);
            weights.resize(nord*nord);

// interpolation matrix is indexed from (1,1) as they are
// automatically generated from Maxima code. Unfortunately, Maxima
// only writes out (1,1) based matrices.
            int start[2] = {1, 1};
            unsigned shape[2];
            shape[0] = nord*nord; shape[1] = nlocal;
            interpMat = Lucee::Matrix<double>(shape, start);
          }

/** Matrix of ordinates */
          Lucee::Matrix<double> ords;
/** Vector of weights */
          std::vector<double> weights;
/** Interpolation matrix */
          Lucee::Matrix<double> interpMat;
      };

/** Data for 2-node Gaussian quadrature */
      GaussQuadData gauss2;
/** Data for 3-node Gaussian quadrature */
      GaussQuadData gauss3;

/**
 * Create matrices for 1st order element.
 */
      void setupPoly1();

/**
 * Create matrices for 2nd order element.
 */
      void setupPoly2();

/**
 * Setup data needed for Gaussian quadrature rules.
 *
 * @param nord Number of ordinates in each direction.
 * @param qData Data to setup.
 */
      void setupGaussQuadData(unsigned nord, GaussQuadData& qData);

/**
 * Helper function to copy data from/to a flat array, given a Lucee
 * field. This method also takes into account the numbering of nodes
 * in the "ghost" cells.
 *
 * @param i I-th index into field.
 * @param j J-th index into field.
 * @param glob On output, these are the exclusive global indices in cell (i,j)
 * @param loc On output, these are the local owned indices (offset 0).
 */
      void getGlobalIndices(int i, int j, std::vector<int>& glob, 
        std::vector<int>& loc);
  };
}

#endif // LC_SERENDEPITY_ELEMENT_2D_H
