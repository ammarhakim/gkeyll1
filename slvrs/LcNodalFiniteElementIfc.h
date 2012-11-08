/**
 * @file	LcNodalFiniteElementIfc.h
 *
 * @brief       Interface class for a reference nodal finite element.
 */

#ifndef LC_NODAL_FINITE_ELEMENT_IFC_H
#define LC_NODAL_FINITE_ELEMENT_IFC_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBasicObj.h>
#include <LcGridIfc.h>
#include <LcMatrix.h>
#include <LcStructGridField.h>

// std includes
#include <cmath>
#include <vector>

namespace Lucee
{
/**
 * Interface class for a reference nodal finite element. This class
 * provides methods (implemented by children) to compute the various
 * things needed to implement the continuous and discontinuous
 * Galerkin finite element methods.
 *
 */
  template <unsigned NDIM>
  class NodalFiniteElementIfc : public Lucee::BasicObj
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
 * Set the current cell location in grid to (i).
 *
 * @param i Index location into grid.
 */
      void setIndex(int i) const;

/**
 * Set the current cell location in grid to (i, j).
 *
 * @param i Index location into grid.
 * @param j Index location into grid.
 */
      void setIndex(int i, int j) const;

/**
 * Set the current cell location in grid to (i, j, k).
 *
 * @param i Index location into grid.
 * @param j Index location into grid.
 * @param k Index location into grid.
 */
      void setIndex(int i, int j, int k) const;

/**
 * Set the current cell location in grid to specified index.
 *
 * @param idx Index location into grid.
 */
      void setIndex(const int idx[NDIM]) const;

/**
 * Get number of local nodes in element.
 *
 * @return number of nodes in element.
 */
      unsigned getNumNodes() const
      {
        return numNodes;
      }

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
 * Method that performs registration of Lua functions.
 *
 * @param lfm Lua function map object.
 */
      static void appendLuaCallableMethods(Lucee::LuaFuncMap& lfm);

/**
 * Lua callable method to get number of exclusive nodes.
 *
 * @param L Lua state to work with.
 * @return number of return values.
 */
      static int luaNumExclusiveNodes(lua_State *L);

/**
 * Lua callable method to get number nodes.
 *
 * @param L Lua state to work with.
 * @return number of return values.
 */
      static int luaNumNodes(lua_State *L);

    protected:
/** Index into current cell */
      mutable int currIdx[NDIM];

/**
 * Create new element with specified number of nodes.
 *
 * @param numNodes Number of nodes in element.
 */
      NodalFiniteElementIfc(unsigned numNodes);

/** 
 * Set number of nodes in element.
 *
 * @param numNodes Number of nodes.
 */
      void setNumNodes(unsigned nN) { numNodes = nN; }

/**
 * Get grid on which elements should be constructured.
 *
 * @return reference to grid.
 */
      template <typename G>
      const G&
      getGrid() const
      {
        if (grid == 0)
          throw Lucee::Except("NodalFiniteElementIfc::getGrid: grid pointer is not valid");
        return dynamic_cast<const G&>(*grid);
      }

    private:
/** Number of nodes */
      unsigned numNodes;
/** Grid on which updater should be applied */
      const Lucee::GridIfc *grid;
  };
}

#endif // LC_NODAL_FINITE_ELEMENT_IFC_H
