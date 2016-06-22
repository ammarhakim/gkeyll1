/**
 * @file	LcPoissonBracketAdvectionEquation1D.h
 *
 * @brief	
 */
#ifndef LC_POISSON_BRACKET_ADVECTION_EQUATION_1D_H
#define LC_POISSON_BRACKET_ADVECTION_EQUATION_1D_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcPoissonBracketEquation.h>
#include <LcAlignedRectCoordSys.h>

namespace Lucee
{
/**
 * Represents the poisson bracket operation for gyrokinetic equation.
 */
  class PoissonBracketAdvectionEquation1D : public Lucee::PoissonBracketEquation
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create a new equation.
 */
      PoissonBracketAdvectionEquation1D();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

      virtual void computeAlphaAtQuadNodes(const Eigen::MatrixXd& hamiltonian, const Eigen::MatrixXd& interpMat,
          const int idx[], Eigen::MatrixXd& alpha);
  };
}

#endif //  LC_POISSON_BRACKET_ADVECTION_EQUATION_1D_H
