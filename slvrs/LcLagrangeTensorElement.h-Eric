/**
 * @file	LcLagrangeTensorElement.h
 *
 * @brief       Lagrange tensor-product element.
 */

#ifndef LC_LAGRANGE_TENSOR_ELEMENT_H
#define LC_LAGRANGE_TENSOR_ELEMENT_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLagrangeTensorBasisCalc.h>
#include <LcNodalFiniteElementIfc.h>

namespace Lucee
{
/**
 * Lagrange tensor-product element. The nodes are laid out in
 * row-major order. Optionally, each direction can have a different
 * polynomial order, and node location can be uniform, Lobatto (both
 * including the element faces) or Gaussian. In the last, there are no
 * nodes on faces, and this should be kept in mind.
 */
  template <unsigned NDIM>
  class LagrangeTensorElement : public Lucee::NodalFiniteElementIfc<NDIM>
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create a new Lagrange tensor element. This does not create a usable
 * object which can only be created from Lua.
 */
      LagrangeTensorElement();

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
 * Extract nodal data at current grid location from field and copy it
 * into a vector. This basically "flattens" the nodal data consistent
 * with the node layout and the stiffness, mass matrices. The output
 * vector should be pre-allocated.
 *
 * @param fld Field to extract nodal data from.
 * @param data On output, this containts a copy of extracted data.
 */
      void extractFromField(const Lucee::Field<NDIM, double>& fld,
        std::vector<double>& data);

/**
 * Copy all nodal data from field and put it into the data array. The
 * data pointer should be pre-allocated.
 *
 * @param fld Field to extract data from.
 * @param data Data space to copy into.
 */
      void copyAllDataFromField(const Lucee::Field<NDIM, double>& fld, double *data);

/**
 * Copy all nodal data to field from a data array.
 *
 * @param data Data space to copy from.
 * @param fld Field to copy data to.
 */
      void copyAllDataToField(const double *data, Lucee::Field<NDIM, double>& fld);

    private:
/** Basis function calculator */
      Lucee::LagrangeTensorBasisCalc<NDIM> basisCalc;
/** Number of 1d volume quadrature nodes to compute */
      unsigned num1DGaussPoints;
/** Number of nodes in each direction */
      unsigned numNodes[NDIM];
/** Indices of exclusively owned nodes */
      std::vector<int> exclusiveNodes;
/** Sequencer for looping over nodes */
      mutable Lucee::RowMajorSequencer<NDIM> nodeSeq;
/** Indexer for global to local mapping */
      Lucee::RowMajorIndexer<NDIM> local2Global;
/** Region for global to local mapping */
      Lucee::Region<NDIM, int> local2GlobalRgn;
/** Strides for use in glocal to local mapping */
      unsigned lgStrides[NDIM];
/** Total number of volume quadrature nodes */
      unsigned numGaussVolNodes;

/**
 * Struct to store list of indexices
 */
      struct IndexList
      {
          std::vector<std::vector<int> > indices;
      };

/** List of node indices on lower surface */
      IndexList lowerIndices[NDIM];
/** List of node indices on upper surface */
      IndexList upperIndices[NDIM];
/** List of exclusively owned node indices */
      IndexList exclusiveNodeIndices;

/** Nodal coordinates relative to lower-left vertex */
      Lucee::Matrix<double> localNodeCoords;
/** Nodal layout */
      Lucee::Node_t nodeLoc;
/** Number of global nodes */
      unsigned numGlobalNodes;
/** Mass matrix */
      Lucee::Matrix<double> mass;
/** Stiffness matrix */
      Lucee::Matrix<double> stiff;
/** Grad-stiffness matrix */
      Lucee::Matrix<double> gradStiff[NDIM];
/** Lower face mass-matrix */
      Lucee::Matrix<double> lowerFace[NDIM];
/** Upper face mass-matrix */
      Lucee::Matrix<double> upperFace[NDIM];

      struct QuadData
      {
/** Interpolation matrix for volume integral */
          Lucee::Matrix<double> interp;
/** Ordinates for volume interpolation */
          Lucee::Matrix<double> ordinates;
/** Weights for volume interpolation */
          std::vector<double> weights;
      };

/** Quadrature data for volume integral */
      QuadData volumeQuad;
/** Quadrature data for surface integral on lower faces */
      QuadData lowerSurfQuad[NDIM];
/** Quadrature data for surface integral on upper faces */
      QuadData upperSurfQuad[NDIM];

/** Weights for nodal quadrature */
      std::vector<double> nodalWeights;
  };
}

#endif // LC_LAGRANGE_TENSOR_ELEMENT_H
