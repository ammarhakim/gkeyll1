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
#include <LcMatrix.h>

namespace Lucee
{
/**
 * Interface class for a reference nodal finite element. This class
 * provides methods (implemented by children) to compute the various
 * things needed to implement the continuous and discontinuous
 * Galerkin finite element methods.
 *
 */
  class NodalFiniteElementIfc : public Lucee::BasicObj
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Get number of nodes in element.
 *
 * @return number of nodes in element.
 */
      unsigned getNumNodes() const
      {
        return numNodes;
      }

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

    private:
/** Number of nodes */
      unsigned numNodes;
  };
}

#endif // LC_NODAL_FINITE_ELEMENT_IFC_H
