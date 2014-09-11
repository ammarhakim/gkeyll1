/**
 * @file	LcPoissonBracketEquation.h
 *
 * @brief	Interface to poisson bracket equations.
 */
#ifndef LC_POISSON_BRACKET_EQUATION_H
#define LC_POISSON_BRACKET_EQUATION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBasicObj.h>
#include <LcConstFieldPtr.h>
#include <LcFieldPtr.h>
#include <LcMatrix.h>
#include <LcRectCoordSys.h>
#include <LcStructGridField.h>

// eigen includes
#include <Eigen/LU>

namespace Lucee
{
  class PoissonBracketEquation : public Lucee::BasicObj
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

      PoissonBracketEquation();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Compute the components of the alpha vector at various quadrature points
 * @param gradHamiltonian gradient of Hamiltonian evaluated at various quadrature points
 *  each row is a different direction, each column is a different quadrature point
 * @param interpMat interpolation matrix for surface or volume. each column is a different
 *  basis function, each row is a different quadrature point
 * @param idx index of aux field
 * @param alpha, each row is a different direction, each column is a different quadrature point
 */
      virtual void computeAlphaAtQuadNodes(const Eigen::MatrixXd& gradHamiltonian, const Eigen::MatrixXd& interpMat,
          const int idx[], Eigen::MatrixXd& alpha);
  };
}

#endif //  LC_POISSON_BRACKET_EQUATION_H
