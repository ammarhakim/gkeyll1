/**
 * @file	LcSerendipityElement.h
 *
 * @brief Serendipity element implemented so far for 2 and 3 dimensions.
 */

#ifndef LC_SERENDIPITY_ELEMENT_H
#define LC_SERENDIPITY_ELEMENT_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcNodalFiniteElementIfc.h>
#include <LcCDIM.h>

// std includes
#include <cmath>
#include <vector>

// eigen includes
#include <Eigen/Dense>
#include <Eigen/LU>

//blitz includes
#include <blitz/array.h>

// etc includes
#include <quadrule.hpp>

namespace Lucee
{
/**
 * Serendipity Elements in arbitrary dimensions
 *
 */
  template <unsigned NDIM>
  class SerendipityElement : public Lucee::NodalFiniteElementIfc<NDIM>
  {
    // Number of components for coordinate arrays etc.
    static const unsigned NC = Lucee::CDIM<NDIM>::N;

    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create a new serendipity element. This does not create a usable
 * object which can only be created from Lua.
 */
      SerendipityElement();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Get local indices of nodes exclusively owned by each cell.
 *
 * @param ndIds On output indices. Vector is cleared and data filled in.
 */
      virtual void getExclusiveNodeIndices(std::vector<int>& ndIds);

/**
 * Get number of surface nodes along lower face in specified
 * direction.
 *
 * @param dir Direction to which face is perpendicular.
 * @return number of nodes in element.
 */
      virtual unsigned getNumSurfLowerNodes(unsigned dir) const;

/**
 * Get number of surface nodes along upper face in specified
 * direction.
 *
 * @param dir Direction to which face is perpendicular.
 * @return number of nodes in element.
 */
      virtual unsigned getNumSurfUpperNodes(unsigned dir) const;

/**
 * Get number of global nodes in element.
 *
 * @return number of nodes in element.
 */
      virtual unsigned getNumGlobalNodes() const;

/**
 * Get mapping of local node numbers in the current cell to global
 * node number. The input vector must be pre-allocated.
 *
 * @param lgMap Local node number to global node number mapping.
 */
      virtual void getLocalToGlobal(std::vector<int>& lgMap) const;

/**
 * Get mapping of local node numbers to global node numbers on lower
 * face of element in direction 'dir' in the current cell. The output
 * vector must be pre-allocated.
 *
 * @param dir Direction to which face is perpendicular.
 * @param lgMap Local node number to global node number mapping.
 */
      virtual void getSurfLowerLocalToGlobal(unsigned dir,
        std::vector<int>& lgMap) const;

/**
 * Get mapping of local node numbers to global node numbers on upper
 * face of element in direction 'dim' in the current cell. The output
 * vector must be pre-allocated.
 *
 * @param lgMap Local node number to global node number mapping.
 */
      virtual void getSurfUpperLocalToGlobal(unsigned dim,
        std::vector<int>& lgMap) const;

/**
 * Get node numbers of the nodes on specified face of element. The
 * output vector must be pre-allocated.
 *
 * @param dir Direction to which face is perpendicular.
 * @param nodeNum Node numbers on face.
 */
      virtual void getSurfLowerNodeNums(unsigned dir,
        std::vector<int>& nodeNum) const;

/**
 * Get node numbers of the nodes on specified face of element. The
 * output vector must be pre-allocated.
 *
 * @param dir Direction to which face is perpendicular.
 * @param nodeNum Node numbers on face.
 */
      virtual void getSurfUpperNodeNums(unsigned dir,
        std::vector<int>& nodeNum) const;

/**
 * Get coordinates of all nodes in element. The output matrix
 * 'nodeCoords' should be pre-allocated have shape numNodes X 3.
 *
 * @param nodeCoords Node coordinates. Should be pre-allocated.
 */
      virtual void getNodalCoordinates(Lucee::Matrix<double>& nodeCoords);

/**
 * Get weights for quadrature. The output vector should be
 * pre-allocated.
 *
 * @param w Weights for quadrature.
 */
      virtual void getWeights(std::vector<double>& w);

/**
 * Get weights for quadrature on upper face. The output vector should
 * be pre-allocated.
 *
 * @param dir Direction to which face is perpendicular.
 * @param w Weights for quadrature.
 */
      virtual void getSurfUpperWeights(unsigned dir,
        std::vector<double>& w);

/**
 * Get weights for quadrature on upper face. The output vector should
 * be pre-allocated.
 *
 * @param dir Direction to which face is perpendicular.
 * @param w Weights for quadrature.
 */
      virtual void getSurfLowerWeights(unsigned dir,
        std::vector<double>& w);

/**
 * Get mass matrix for this reference element. The output matrix
 * should be pre-allocated.
 *
 * @param NjNk On output, mass matrix of element.
 */
      virtual void getMassMatrix(Lucee::Matrix<double>& NjNk) const;

/**
 * Get mass matrix for this reference element. The output matrix
 * should be pre-allocated.
 *
 * @param dir Direction of face.
 * @param NjNk On output, mass matrix of element.
 */
      virtual void getLowerFaceMassMatrix(unsigned dir,
        Lucee::Matrix<double>& NjNk) const;

/**
 * Get mass matrix for this reference element. The output matrix
 * should be pre-allocated.
 *
 * @param dir Direction of face.
 * @param NjNk On output, mass matrix of element.
 */
      virtual void getUpperFaceMassMatrix(unsigned dir,
        Lucee::Matrix<double>& NjNk) const;

/**
 * Get stiffness matrix (grad.Nj \dot grad.Nk) for this reference
 * element. The output matrix should be pre-allocated.
 *
 * @param DNjDNk On output, stiffness matrix of element.
 */
      virtual void getStiffnessMatrix(Lucee::Matrix<double>& DNjDNk) const;

/**
 * Get partial stiffness matrix (grad.Nj Nk) for this reference
 * element. The output matrix should be pre-allocated.
 *
 * @param dir Direction for gradient.
 * @param DNjNk On output, partial stiffness matrix of element.
 */
      virtual void getGradStiffnessMatrix(unsigned dir, Lucee::Matrix<double>& DNjNk) const;

/**
 * Get number of nodes needed for Gaussian quadrature in the element
 * interior.
 *
 * @return Number of nodes needed for Gaussian quadrature.
 */
      virtual unsigned getNumGaussNodes() const;

/**
 * Get number of nodes needed for Gaussian quadrature on the element
 * surface.
 *
 * @return Number of nodes needed for Gaussian quadrature.
 */
      virtual unsigned getNumSurfGaussNodes() const;

/**
 * Get data needed for Gaussian quadrature for this element. All
 * output matrices and vectors must be pre-allocated.
 *
 * @param interpMat On output, interpolation matrix.
 * @param ordinates On output, quadrature ordinates.
 * @param weights On output, quadrature weights.
 */
      virtual void getGaussQuadData(Lucee::Matrix<double>& interpMat,
        Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const;

/**
 * Get data needed for Gaussian quadrature on lower surfaces of this
 * element. All output matrices and vectors must be pre-allocated.
 *
 * @param interpMat On output, interpolation matrix.
 * @param ordinates On output, quadrature ordinates.
 * @param weights On output, quadrature weights.
 */
      virtual void getSurfLowerGaussQuadData(unsigned dir, Lucee::Matrix<double>& interpMat,
        Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const;

/**
 * Get data needed for Gaussian quadrature on upper surfaces of this
 * element. All output matrices and vectors must be pre-allocated.
 *
 * @param interpMat On output, interpolation matrix.
 * @param ordinates On output, quadrature ordinates.
 * @param weights On output, quadrature weights.
 */
      virtual void getSurfUpperGaussQuadData(unsigned dir, Lucee::Matrix<double>& interpMat,
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
      virtual void getMomentMatrix(unsigned p, Lucee::Matrix<double>& momMat) const;
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
 * Get matrices needed to compute hyper-diffusion operator. The
 * matrices are for the current cell and each of its face neighbors,
 * stored in "lowerMat" for cells sharing lower faces and "upperMat"
 * for cells sharing upper faces. A linear combination of these
 * matrices when multiplied by the nodal data in the corresponding
 * cells should give the discrete diffusion operator.
 *
 * @param iMat Matrix for current cell.
 * @param lowerMat Matrices for cells sharing lower faces.
 * @param upperMat Matrices for cells sharing upper faces.
 */
      void getHyperDiffusionMatrices(std::vector<Lucee::Matrix<double> >& iMat, 
        std::vector<Lucee::Matrix<double> >& lowerMat, std::vector<Lucee::Matrix<double> >& upperMat) const;

/**
 * Get coefficients for applying reflecting boundary conditions on
 * lower side in direction 'dir'. The vector nodeMap[i] is the
 * reflected node corresponding to node 'i'.
 *
 * @param dir Direction to which face is perpendicular.
 * @param nodeMap Map for reflecting nodes.
 */
      virtual void getLowerReflectingBcMapping(unsigned dir,
        std::vector<unsigned>& nodeMap) const;

/**
 * Get coefficients for applying reflecting boundary conditions on
 * upper side in direction 'dir'. The vector nodeMap[i] is the
 * reflected node corresponding to node 'i'.
 *
 * @param dir Direction to which face is perpendicular.
 * @param nodeMap Map for reflecting nodes.
 */
      virtual void getUpperReflectingBcMapping(unsigned dir,
        std::vector<unsigned>& nodeMap) const;
/**
 * Compute mapping of nodes from a face to element interior for lower
 * face in specified direction. The shape of the output matrix
 * faceToIntMap is (numVolNodes X numFaceNodes).
 *
 * @param dir Direction to which face is perpendicular.
 * @param faceToIntMap Matrix storing the mapping.
 */
      virtual void getLowerFaceToInteriorMapping(unsigned dir,
        Lucee::Matrix<double>& faceToIntMap) const;
/**
 * Extract nodal data at current grid location from field and copy it
 * into a vector. This basically "flattens" the nodal data consistent
 * with the node layout and the stiffness, mass matrices. The output
 * vector should be pre-allocated.
 *
 * @param fld Field to extract nodal data from.
 * @param data On output, this containts a copy of extracted data.
 */
      virtual void extractFromField(const Lucee::Field<NDIM, double>& fld,
        std::vector<double>& data);

/**
 * Copy all nodal data from field and put it into the data array. The
 * data pointer should be pre-allocated.
 *
 * @param fld Field to extract data from.
 * @param data Data space to copy into.
 */
      virtual void copyAllDataFromField(const Lucee::Field<NDIM, double>& fld,
        double *data);

/**
 * Copy all nodal data to field from a data array.
 *
 * @param data Data space to copy from.
 * @param fld Field to copy data to.
 */
      virtual void copyAllDataToField(const double *data, Lucee::Field<NDIM, double>& fld);

/**
 * Evaluate basis functions at location. The results should be stored
 * in the pre-allocated 'vals' vector.
 *
 * @param xc Coordinates in element.
 * @param vals Values of basis functions. Pre-allocated.
 */
      virtual void evalBasis(double xc[NDIM], std::vector<double>& vals) const;
  
  private :
/** Polynomial order of element */
      unsigned polyOrder;
/** Maximum polynomial power anticipated */
      unsigned maxPower;
/** Total number of quadrature points in 1-D */
      unsigned numGaussPoints;
/** Highest moment degree to compute */
      int maxMoment;
/** Grid spacing in various dimensions */
      double dq[NDIM];
/** Grid spacing squared in various dimensions */
      double dq2[NDIM];
/** Matrix containing coordinates of node on reference element. Rows = nodes, Cols = dim */
      Eigen::MatrixXd nodeList;
      std::vector<blitz::Array<double,NDIM> > functionVector;
/** Matrix containing basis functions evaluated at volume gaussian integration locations
    Correspondance between column and gaussian node set is kept track of in gaussNodeList
    Each row is a different quadrature point.
    Each column is a different basis function evaluated at the same point*/
      Eigen::MatrixXd functionEvaluations;
/** Matrix containing basis functions evaluated at surface gaussian integration locations
    Correspondance between row and gaussian node set is kept track of in gaussNodeListUpperSurf
    Each row is a different evaluation location
    Each column is a different basis function */
      std::vector<Eigen::MatrixXd> upperSurfaceEvaluations;
/** Matrix containing basis functions evaluated at surface gaussian integration locations
    Correspondance between row and gaussian node set is kept track of in gaussNodeListLowerSurf
    Each row is a different evaluation location
    Each column is a different basis function */
      std::vector<Eigen::MatrixXd> lowerSurfaceEvaluations;
/** Matrix whose ith row corresponds to the coordinates of the node at which the ith column of functionEvaluations
    (and functionDEvaluations) is evaluated at. Size NDIM+1 columns because last entry is net weight */
      Eigen::MatrixXd gaussNodeList;
/** Vector of gaussian quadrature coordinates on upper surface indexed by dimension
    Each row is the coordinates of a node + weight in last col */
      std::vector<Eigen::MatrixXd> gaussNodeListUpperSurf;
/** Vector of gaussian quadrature coordinates on lower surface indexed by dimension*/
      std::vector<Eigen::MatrixXd> gaussNodeListLowerSurf;
/** Weights for quadrature (one dimension)*/
      std::vector<double> gaussWeights;
/** Ordinates for (one dimension) quadrature */
      std::vector<double> gaussPoints;
/** Vector of face-mass matrices indexed by dimension (Lower) */
      std::vector<Eigen::MatrixXd> refFaceMassLower;
/** Vector of face-mass matrices indexed by dimension (Upper) */
      std::vector<Eigen::MatrixXd> refFaceMassUpper;
/** Vector of grad stiffness matrices indexed by dimension */
      std::vector<Eigen::MatrixXd> refGradStiffness;
/** Mass matrix in reference coordinates */
      Eigen::MatrixXd refMass;
/** Stiffness matrix in reference coordinates */
      Eigen::MatrixXd refStiffness;
/** Vector of moment matrices indexed by moment value p  */
      std::vector<Eigen::MatrixXd> momMatrix;
/** List of matrices for current cell */
      std::vector<Eigen::MatrixXd> iMatDiffusion;
      std::vector<Eigen::MatrixXd> iMatHyperDiffusion;
/** List of matrices on each lower face */
      std::vector<Eigen::MatrixXd> lowerMatDiffusion;
      std::vector<Eigen::MatrixXd> lowerMatHyperDiffusion;
/** List of matrices on each upper face */
      std::vector<Eigen::MatrixXd> upperMatDiffusion;
      std::vector<Eigen::MatrixXd> upperMatHyperDiffusion;
/** Face to interior mapping matrices */
      std::vector<Eigen::MatrixXd> lowerFaceToInteriorMapMatrices;
/**
 *    Create necessary matrices needed for 1,2,3rd order serendipity elements.
 *    Currently only works for 3-D cases.
 */
      void setupMatrices();
/**
 *    Resize output matrices computed in setupMatrices()
 */
      void resizeMatrices();
/**
 *    Populate nodeList with serendipity node locations on reference
 *    element
 */
      void getNodeList(Eigen::MatrixXd& nodeMatrix, int degree);
/**
 *    Create basis monomials by populating matrix of rows
 *    [a b c] to represent x^a*y^b*z^c monomials
 */
      void setupBasisMatrix(Eigen::MatrixXi& basisMatrix, int degree);
/**
  *   Compute the basis functions in terms of the basis monomials
  */
      void computeBasisFunctions(std::vector<blitz::Array<double, NDIM> >& functionVector,
        const Eigen::MatrixXd& nodeList, int degree);
 /**
 *    Evaluate a polynomial represented by coefficients in a n-d array at a specific location
 *    defined by a vector nodeCoords
 */
      double evalPolynomial(const blitz::Array<double, NDIM>& polyCoeffs, const Eigen::VectorXd& nodeCoords) const;

/**
*     Compute the partial derivative of a polynomial in direction 'dir'
*/
      blitz::Array<double, NDIM> computePolynomialDerivative(const blitz::Array<double, NDIM>& poly, int dir);

/**
*     Compute the mass matrix on the reference element.
*/
      void computeMass(Eigen::MatrixXd& resultMatrix);

/**
*     Compute the face-mass matrices on the reference element in direction num.
*/
      void computeFaceMass(int dir, Eigen::MatrixXd& lowerResultMatrix, Eigen::MatrixXd& upperResultMatrix);

/**
*     Compute the stiffness matrix on the reference element.
*/
      void computeStiffness(const blitz::Array<double, 3>& functionDerivative, Eigen::MatrixXd& resultMatrix);     

/**
 *     Compute the grad stiffness matrix in direction dir on the reference element.
 */
      void computeGradStiffness(const blitz::Array<double, 3>& functionDerivative,
        int dir, Eigen::MatrixXd& resultMatrix);
/**
 *    Compute 'cross-dimensional' basis functions
 */
      void setupMomentMatrices();
/**
 *     Compute the det(J) factor
 */
      void computeTransformationScales();
/**
 *    Compute the number of degrees of freedom for serendipity element
 *    with degree r and in dimension n
 */
      int getSerendipityDimension(int degree, int dimension) const;
/**
 *    Compute n!
 */
      int factorial(int n) const;
  };
}

#endif // LC_SERENDIPITY_ELEMENT_H
