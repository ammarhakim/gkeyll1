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
 *             3          2
 *             o----------o
 *             |          |
 *             |          |
 *             |          |
 *             |          |
 *             o----------o
 *             0          1
 *
 * Order 2 element
 *
 *             3     6     2
 *             o-----o-----o
 *             |           |
 *             |           |
 *           7 o           o 5
 *             |           |
 *             |           |
 *             o-----o-----o
 *             0     4     1
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
      void getExclusiveNodeIndices(std::vector<int>& ndIds);

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
 * Get number of nodes needed for Gaussian quadrature in the element
 * interior.
 *
 * @return Number of nodes needed for Gaussian quadrature.
 */
      unsigned getNumGaussNodes() const;

/**
 * Get number of nodes needed for Gaussian quadrature on the element
 * surface.
 *
 * @return Number of nodes needed for Gaussian quadrature.
 */
      unsigned getNumSurfGaussNodes() const;

/**
 * Get data needed for Gaussian quadrature of specified order on this
 * element. All output matrices and vectors must be pre-allocated.
 *
 * @param interpMat On output, interpolation matrix.
 * @param ordinates On output, quadrature ordinates (numNodes X 3)
 * @param weights On output, quadrature weights.
 */
      void getGaussQuadData(Lucee::Matrix<double>& interpMat,
        Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const;

/**
 * Get data needed for Gaussian quadrature on lower surfaces of this
 * element. All output matrices and vectors must be pre-allocated.
 *
 * @param interpMat On output, interpolation matrix.
 * @param ordinates On output, quadrature ordinates.
 * @param weights On output, quadrature weights.
 */
      void getSurfLowerGaussQuadData(unsigned dir, Lucee::Matrix<double>& interpMat,
        Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const;

/**
 * Get data needed for Gaussian quadrature on upper surfaces of this
 * element. All output matrices and vectors must be pre-allocated.
 *
 * @param interpMat On output, interpolation matrix.
 * @param ordinates On output, quadrature ordinates.
 * @param weights On output, quadrature weights.
 */
      void getSurfUpperGaussQuadData(unsigned dir, Lucee::Matrix<double>& interpMat,
        Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const;

/**
 * Get matrix for projection on a lower-dimensional basis set. This
 * method returns the moment matrix (in 2D, for example)
 *
 * \int y^p \phi(x,y) \psi(x) dx dy
 *
 * @param p Required moment.
 * @param momMatt On output, moment matrix.
 */
      void getMomentMatrix(unsigned p, Lucee::Matrix<double>& momMat) const;

/**
 * Get matrices needed to compute diffusion operator. The matrices are
 * for the current cell and each of its face neighbors, stored in
 * "lowerMat" for cells sharing lower faces and "upperMat" for cells
 * sharing upper faces. A linear combination of these matrices when
 * multiplied by the nodal data in the corresponding cells should give
 * the discrete diffusion operator.
 *
 * @param iMat Matrix for current cell, split into contributions from each direction.
 * @param lowerMat Matrices for cells sharing lower faces.
 * @param upperMat Matrices for cells sharing upper faces.
 */
      void getDiffusionMatrices(std::vector<Lucee::Matrix<double> >& iMat,
        std::vector<Lucee::Matrix<double> >& lowerMat, std::vector<Lucee::Matrix<double> >& upperMat) const;

/**
 * Get coefficients for applying reflecting boundary conditions on
 * lower side in direction 'dir'. The vector nodeMap[i] is the
 * reflected node corresponding to node 'i'.
 *
 * @param dir Direction to which face is perpendicular.
 * @param nodeMap Map for reflecting nodes.
 */
      void getLowerReflectingBcMapping(unsigned dir, std::vector<unsigned>& nodeMap) const;

/**
 * Get coefficients for applying reflecting boundary conditions on
 * upper side in direction 'dir'. The vector nodeMap[i] is the
 * reflected node corresponding to node 'i'.
 *
 * @param dir Direction to which face is perpendicular.
 * @param nodeMap Map for reflecting nodes.
 */
      void getUpperReflectingBcMapping(unsigned dir, std::vector<unsigned>& nodeMap) const;

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
/** List of diffusion matrices for current cell */
      std::vector<Lucee::Matrix<double> > iMatDiffusion;
/** List of diffusion matrices on each lower face */
      std::vector<Lucee::Matrix<double> > lowerMatDiffusion;
/** List of diffusion matrices on each upper face */
      std::vector<Lucee::Matrix<double> > upperMatDiffusion;

/** Mapping for reflection on lower faces in X */
      std::vector<unsigned> lowerNodeMap0;
/** Mapping for reflection on lower faces in Y */
      std::vector<unsigned> lowerNodeMap1;
/** Mapping for reflection on upper faces in X */
      std::vector<unsigned> upperNodeMap0;
/** Mapping for reflection on upper faces in Y */
      std::vector<unsigned> upperNodeMap1;

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
          GaussQuadData(unsigned nord = 1, unsigned nlocal = 1)
            : ords(nord, 3), weights(nord), interpMat(nord, nlocal)
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
            ords = Lucee::Matrix<double>(nord, 3);
            weights.resize(nord);

// interpolation matrix is indexed from (1,1) as they are
// automatically generated from Maxima code. Unfortunately, Maxima
// only writes out (1,1) based matrices.
            int start[2] = {1, 1};
            unsigned shape[2];
            shape[0] = nord; shape[1] = nlocal;
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
/** Quadrature, lower faces */
      GaussQuadData gauss2Lower[2];
/** Quadrature, upper faces */
      GaussQuadData gauss2Upper[2];

/** Data for 3-node Gaussian quadrature */
      GaussQuadData gauss3;
/** Quadrature, lower faces */
      GaussQuadData gauss3Lower[2];
/** Quadrature, upper faces */
      GaussQuadData gauss3Upper[2];

/**
 * Struct to hold moment matrices.
 */
      struct MomentMatrix
      {
/**
 * Create moment matrix.
 */
          MomentMatrix()
            : m(2,4)
          {
          }

/** Moment matrix */
          Lucee::Matrix<double> m;
      };

/** List of moment matrices */
      MomentMatrix momMatrix[3]; // for now 0, 1 and 2nd moments are computed

/**
 * Create matrices for 1st order element.
 */
      void setupPoly1();

/**
 * Create matrices for 2nd order element.
 */
      void setupPoly2();

/**
 * Setup data needed for Gaussian quadrature rules for polynomial
 * order 1.
 */
      void setupGaussQuadDataPoly1();

/**
 * Setup data needed for Gaussian quadrature rules for polynomial
 * order 2.
 */
      void setupGaussQuadDataPoly2();

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
