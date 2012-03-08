/**
 * @file	LcLobattoElement1D.h
 *
 * @brief       Reference finite element with Lobatto nodes
 */

#ifndef LC_LOBATTO_ELEMENT_1D_H
#define LC_LOBATTO_ELEMENT_1D_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcNodalFiniteElementIfc.h>

namespace Lucee
{
/**
 * Lobatto element in 2D. The reference element is the interval [-1,1]
 * X [-1,1].
 */
  class LobattoElement1D : public Lucee::NodalFiniteElementIfc<1>
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create a new serendipity element. This does not create a usable
 * object which can only be created from Lua.
 */
      LobattoElement1D();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Get number of global nodes in element.
 *
 * @return number of nodes in element.
 */
      unsigned getNumGlobalNodes() const;

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
      void getMassMatrix(Lucee::Matrix<double> NjNk) const;

/**
 * Get stiffness matrix (grad.Nj \dot grad.Nk) for this reference
 * element. The output matrix should be pre-allocated.
 *
 * @param DNjDNk On output, stiffness matrix of element.
 */
      void getStiffnessMatrix(Lucee::Matrix<double> DNjDNk) const;

    private:
/** Order of polynomial in element */
      unsigned polyOrder;
/** Mass matrix in reference coordinates */
      Lucee::Matrix<double> refNjNk;
/** Stiffness matrix in reference coordinates */
      Lucee::Matrix<double> refDNjDNk;
/** Total global degrees of freedom */
      unsigned numGlobal;

/**
 * Create matrices for polyOrder 1.
 */
      void setupPoly1();

/**
 * Create matrices for polyOrder 2.
 */
      void setupPoly2();

/**
 * Create matrices for polyOrder 3.
 */
      void setupPoly3();
  };
}

#endif // LC_LOBATTO_ELEMENT_1D_H
