/**
 * @file	LcLagrangeTensorBasisCalc.h
 *
 * @brief	Class to calculate data needed for Lagrange tensor basis functions.
 */
#ifndef LC_LAGRANGE_TENSOR_BASIS_CALC_H
#define LC_LAGRANGE_TENSOR_BASIS_CALC_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// gkeyll includes
#include <LcMatrix.h>
#include <LcRowMajorIndexer.h>

// std includes
#include <vector>

namespace Lucee
{
/** Enum for location of nodes */
  enum Node_t { LOBATTO, GAUSSIAN, UNIFORM };

/**
 * Class to calculate data needed for Lagrange tensor basis functions.
 */
  template <unsigned NDIM>
  class LagrangeTensorBasisCalc
  {
    public:
/** We need to friend ourself to allow accessing private stuff from another dimension */
      template <unsigned RDIM> friend class LagrangeTensorBasisCalc;

/**
 * Create new tensor basis calculator object.
 */
      LagrangeTensorBasisCalc();

/**
 * Get total number of nodes in element.
 *
 * @return Total number of nodes in element.
 */
      unsigned getNumNodes() const 
      { return totalNodes; }

/**
 *  Calculate the data needed in constructing the basis
 *  function. Calling this function again will redo the calculation.
 *
 * @param type Type of nodal layout.
 * @param numNodes Number of nodes in each direction. Should have exactly ndim elements.
 */
      void calc(Node_t type, const unsigned numNodes[NDIM]);

/**
 * Get node location in specified direction.
 *
 * @param dir Direction for node locations.
 * @return vector of node locations.
 */
      std::vector<double> getNodeLoc(unsigned dir) const
      {
        return nodeLocs[dir].loc;
      }

/**
 * Get coordinate of specified node.
 *
 * @param nIdx Node index (0-based)
 * @param xn On output, coordinate of node.
 */
      void fillWithNodeCoordinate(unsigned nIdx, double xn[NDIM]) const
      {
        for (unsigned d=0; d<NDIM; ++d)
          xn[d] = nodeCoords[nIdx].x[d];
      }

/**
 * Fetch coefficient matrix.
 *
 * @param coeff On output, this contains the coefficient matrix. Should be pre-allocated.
 */
      void getCoeffMat(Lucee::Matrix<double>& coeff) const
      {
        coeff.copy(expandCoeff);
      }

/**
 * Fetch mass-matrix.
 *
 * @param mMatrix On output, this contains the mass matrix. Should be pre-allocated.
 */
      void getMassMatrix(Lucee::Matrix<double>& mMatrix) const
      {
        mMatrix.copy(massMatrix);
      }

/**
 * Fetch stiffness-matrix.
 *
 * @param dx Grid spacing in each direction
 * @param mMatrix On output, this contains the stiffness matrix. Should be pre-allocated.
 */
      void getStiffnessMatrix(double dx[NDIM], Lucee::Matrix<double>& mMatrix) const
      {
        for (unsigned r=0; r<totalNodes; ++r)
        {
          for (unsigned c=0; c<totalNodes; ++c)
          {
            mMatrix(r,c) = 0.0;
            for (unsigned d=0; d<NDIM; ++d)
              mMatrix(r,c) += 4.0/(dx[d]*dx[d])*stiffMatrix[d](r,c);
          }
        }
      }

/**
 * Fetch grad-stiffness matrix in specified direction.
 *
 * @param dir Direction in which matrix is required.
 * @param gMatrix On output, this contains the grad-stiff matrix. Should be pre-allocated.
 */
      void getGradStiffnessMatrix(unsigned dir, Lucee::Matrix<double>& gMatrix) const
      {
        gMatrix.copy(gradStiff[dir]);
      }

/**
 * Get face-mass matrix on lower face in specified direction for element.
 *
 * @param dir Direction of face.
 * @param fMatrix On output, face-mass matrix of element. Should be pre-allocated
 */
      void getLowerFaceMassMatrix(unsigned dir, Lucee::Matrix<double>& fMatrix) const
      {
        fMatrix.copy(lowerFaceMass[dir]);
      }

/**
 * Get face-mass matrix on upper face in specified direction for element.
 *
 * @param dir Direction of face.
 * @param fMatrix On output, face-mass matrix of element. Should be pre-allocated
 */
      void getUpperFaceMassMatrix(unsigned dir, Lucee::Matrix<double>& fMatrix) const
      {
        fMatrix.copy(upperFaceMass[dir]);
      }

/**
 * Get data needed for Gaussian quadrature for this element. All
 * output matrices and vectors must be pre-allocated.
 *
 * @param interpMat On output, interpolation matrix.
 * @param ordinates On output, quadrature ordinates.
 * @param weights On output, quadrature weights.
 */
      void getGaussQuadData(Lucee::Matrix<double>& interpMat,
        Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const
      {
        interpMat.copy(volumeGaussInterp);
        ordinates.copy(volumeGaussOrdinates);
        weights.clear(); weights.resize(totalNodes);
        weights = volumeGaussWeights;
      }

/**
 * Evaluate specified basis function at location. This method is very
 * slow, and so should not be called inside an inner loop. It is
 * provided to allow initialization of data needed in updaters, so its
 * cost will be amortized.
 *
 * @param bIdx Basis function index (0-based).
 * @param xc Coorinates in element.
 * @return value of basis function at location.
 */
      double evalBasis(unsigned bIdx, double xc[NDIM]) const;

/**
 * Get indexer for use in stepping over nodes in the proper order.
 *
 * @return Indexer for stepping over nodes.
 */
      Lucee::RowMajorIndexer<NDIM> getIndexer() const;

/**
 * Get number of surface nodes along lower face in specified
 * direction.
 *
 * @param dir Direction to which face is perpendicular.
 * @return number of nodes in element.
 */
      unsigned getNumSurfLowerNodes(unsigned dir) const
      { 
        return lowerNodes[dir].nodes.size();
      }

/**
 * Get number of surface nodes along upper face in specified
 * direction.
 *
 * @param dir Direction to which face is perpendicular.
 * @return number of nodes in element.
 */
      unsigned getNumSurfUpperNodes(unsigned dir) const
      {
        return upperNodes[dir].nodes.size();
      }

/**
 * Get local indices of nodes exclusively owned by each cell.
 *
 * @param ndIds On output indices. Vector is cleared and data filled in.
 */
      void getExclusiveNodeIndices(std::vector<int>& ndIds)
      {
        ndIds = exclNodes;
      }

/**
 * Get NDIM indices of nodes exclusively owned by each cell.
 *
 * @param ndIds On output indices. Vector is cleared and data filled in.
 */
      void getExclusiveNodeNdimIndices(std::vector<std::vector<int> >& ndIds)
      {
        ndIds.clear();
        ndIds.resize(exclNodesIndices.size());
        for (unsigned i=0; i<exclNodesIndices.size(); ++i)
          ndIds[i] = exclNodesIndices[i];
      }

/**
 * Get node numbers of the nodes on specified face of element. The
 * output vector must be pre-allocated.
 *
 * @param dir Direction to which face is perpendicular.
 * @param nodeNum Node numbers on face.
 */
      void getSurfLowerNodeNums(unsigned dir, std::vector<int>& nodeNum) const
      {
        nodeNum = lowerNodes[dir].nodes;
      }

/**
 * Get node numbers of the nodes on specified face of element. The
 * output vector must be pre-allocated.
 *
 * @param dir Direction to which face is perpendicular.
 * @param nodeNum Node numbers on face.
 */
      void getSurfUpperNodeNums(unsigned dir, std::vector<int>& nodeNum) const
      {
        nodeNum = upperNodes[dir].nodes;
      }

    private:
/** Node layout */
      Node_t nodeLayout;
/** Number of nodes in each direction */
      unsigned numNodes[NDIM];
/** Total number of nodes */
      unsigned totalNodes;
/** Region representing node shape */
      typename Lucee::Region<NDIM, int> nodeRgn;

/** Structure to store location of nodes */
      struct NodeLoc
      {
/** Location of nodes */
          std::vector<double> loc;
      };

/** Location of nodes in each direction */
      NodeLoc nodeLocs[NDIM];

/** Structure to store coordinates of each node */
      struct NodeCoord
      {
/** Node coordinate */
          double x[NDIM];
      };

/** Coordinates of each node (this seems redundant, but is stored
 * anyway to make life easier) */
      std::vector<NodeCoord> nodeCoords;

/** Matrix of expansion coefficients */
      Lucee::Matrix<double> expandCoeff;
/** Mass matrix for element */
      Lucee::Matrix<double> massMatrix;
/** Grad-stiffness matrix */
      Lucee::Matrix<double> gradStiff[NDIM];
/** Stiffness matrix for element: strangely, this is stored
 * per-direction to allow multiplication by grid space when computing
 * the full stiffness matrix */
      Lucee::Matrix<double> stiffMatrix[NDIM];

/** Number of exclusively owned nodes */
      unsigned numExclNodes;
/** Indices of exclusively owned nodes */
      std::vector<int> exclNodes;
/** NDIM Indices of exclusively owned nodes */
      std::vector<std::vector<int> > exclNodesIndices;

/** Structure to hold node list for surface nodes */
      struct NodeList
      {
/** List of nodes */
          std::vector<int> nodes;
      };

/** List of nodes on lower faces */
      NodeList lowerNodes[NDIM];
/** List of nodes on upper faces */
      NodeList upperNodes[NDIM];

/** Lower face mass matrices */
      Lucee::Matrix<double> lowerFaceMass[NDIM];
/** Upper face mass matrices */
      Lucee::Matrix<double> upperFaceMass[NDIM];

/** Interpolation matrix for volume integral */
      Lucee::Matrix<double> volumeGaussInterp;
/** Ordinates for volume interpolation */
      Lucee::Matrix<double> volumeGaussOrdinates;
/** Weights for volume interpolation */
      std::vector<double> volumeGaussWeights;

/**
 * Create nodes located at Lobatto quadrature points.
 */
      void createLobattoNodes();

/**
 * Create nodes located at Gaussian quadrature points.
 */
      void createGaussianNodes();

/**
 * Create nodes with uniform spacing.
 */
      void createUniformNodes();

/**
 * Compute mass matrix.
 */
      void calcMassMatrix();

/**
 * Compute grad-stiffness matrix in specified direction.
 *
 * @param dir Direction to compute grad-stiffness matrix in.
 */
      void calcGradStiff(unsigned dir);

/**
 * Compute lower and upper face mass-matrices in specified direction.
 *
 * @param dir Direction to compute face mass-matrix in.
 */
      void calcFaceMass(unsigned dir);

/**
 * Compute stiffness matrix.
 */
      void calcStiffMatrix();

/**
 * Compute nodal layout.
 */
      void calcBasicData(Node_t type, const unsigned nn[NDIM]);

/**
 * Compute matrix int_{-1}^{1} P_n'(x) P_m'(x) dx.
 *
 * @param dpdp On output, this contains the needed matrix.
 */
      void calcDpDp(Lucee::Matrix<double>& dpdp);

/**
 * Compute volume quadrature data.
 */
      void calcVolumeQuad();
  };

// Explicitly instantiate NDIM=0 case to void compiler barfs. This is
// a rather horrible hack but I do not at present see a way around
// this. (Ammar Hakim, November 19th 2012)
  template <>
  class LagrangeTensorBasisCalc<0>
  {
    public:
      void calcFaceMass(unsigned dir)
      {
// deliberately empty
      }

      void calcBasicData(Node_t type, const unsigned nn[])
      {
// deliberately empty
      }

      void calcMassMatrix()
      {
// deliberately empty
      }

      unsigned getNumNodes() const
      { return 0; }

      void getMassMatrix(Lucee::Matrix<double>& mMatrix) const
      {
// deliberately empty
      }
  };
}

#endif // LC_LAGRANGE_TENSOR_BASIS_CALC_H
