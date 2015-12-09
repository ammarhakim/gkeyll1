/**
 * @file	LcPoissonBracketVlasovPoissonEquation3D.h
 *
 * @brief	Poisson bracket for Vlasov-Poisson problem with 1D2V and static constant B.
 */
#ifndef LC_POISSON_BRACKET_VLASOV_POISSON_3D_H
#define LC_POISSON_BRACKET_VLASOV_POISSON_3D_H

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
  class PoissonBracketVlasovPoissonEquation3D : public Lucee::PoissonBracketEquation
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create a new equation.
 */
      PoissonBracketVlasovPoissonEquation3D();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

      virtual void computeAlphaAtQuadNodes(const Eigen::MatrixXd& hamiltonian, const Eigen::MatrixXd& interpMat,
          const int idx[], Eigen::MatrixXd& alpha);
    private:
      // Mass of species we are updating
      double speciesMass;
      // Charge of species we are updating
      double speciesCharge;
      // Value of background constant and static magnetic field
      double BZ0;
  };
}

#endif //  LC_POISSON_BRACKET_VLASOV_POISSON_3D_H
