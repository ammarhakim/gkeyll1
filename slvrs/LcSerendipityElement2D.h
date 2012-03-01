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

namespace Lucee
{
/**
 * Serendipity reference element in 2D. The reference element is a
 * square [-1,1] X [-1,1]. The supported elements are shown below with
 * the corresponding node layout.
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
  class SerendipityElement2D : public Lucee::NodalFiniteElementIfc
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

/**
 * Create matrices for 1st order element.
 */
      void setupOrder1();

/**
 * Create matrices for 2nd order element.
 */
      void setupOrder2();
  };
}

#endif // LC_SERENDEPITY_ELEMENT_2D_H
