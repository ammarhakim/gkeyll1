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
 * Get number of global nodes in element.
 *
 * @return number of nodes in element.
 */
      virtual unsigned getNumGlobalNodes() const;

/**
 * Evaluate 'n'th basis function at location (x,y) in the reference
 * element.
 *
 * @param n Basis at node 'n'
 * @param x X-location in reference element.
 * @param y Y-location in reference element.
 * @return Value of basis function 'n' at (x,y)
 */
      virtual double evalBasis(unsigned n, double x, double y) const;

/**
 * Get mapping of local node numbers in the current cell to global
 * node number. The input vector must be pre-allocated.
 *
 * @param lgMap Local node number to global node number mapping.
 */
      virtual void getLocalToGlobal(std::vector<int>& lgMap) const;

/**
 * Get mass matrix for this reference element. The output matrix
 * should be pre-allocated.
 *
 * @param NjNk On output, mass matrix of element.
 */
      virtual void getMassMatrix(Lucee::Matrix<double> NjNk) const;

/**
 * Get stiffness matrix (grad.Nj \dot grad.Nk) for this reference
 * element. The output matrix should be pre-allocated.
 *
 * @param DNjDNk On output, stiffness matrix of element.
 */
      virtual void getStiffnessMatrix(Lucee::Matrix<double> DNjDNk) const;

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
