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
 * Serenditpity reference element in 2D.
 */
  class SerendipityElement2D : public Lucee::NodalFiniteElementIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create a new serendipity element.
 *
 * @param order of element.
 */
      SerendipityElement2D(unsigned order);

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
/** Order of element */
      unsigned elmOrder;
/** Mass matrix in reference coordinates */
      Lucee::Matrix<double> refNjNk;
/** Stiffness matrix in reference coordinates */
      Lucee::Matrix<double> refDNjDNk;

/**
 * Create matrices for 2nd order element.
 */
      void setupOrder2();

/**
 * Create matrices for 2nd order element.
 */
      void setupOrder3();
  };
}

#endif // LC_SERENDEPITY_ELEMENT_2D_H
