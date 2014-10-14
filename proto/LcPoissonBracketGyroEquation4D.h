/**
 * @file	LcPoissonBracketGyroEquation4D.h
 *
 * @brief	Compute Jacobian*alpha at all quadrature points
 */
#ifndef LC_POISSON_BRACKET_GYRO_EQUATION_4D_H
#define LC_POISSON_BRACKET_GYRO_EQUATION_4D_H

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
  class PoissonBracketGyroEquation4D : public Lucee::PoissonBracketEquation
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create a new equation.
 */
      PoissonBracketGyroEquation4D();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

      virtual void computeAlphaAtQuadNodes(const Eigen::MatrixXd& hamiltonian, const Eigen::MatrixXd& interpMat,
          const int idx[], Eigen::MatrixXd& alpha);
    private:
      Lucee::Field<4, double> *bStarYField;
      // Mass of species we are updating
      double speciesMass;
      // Charge of species we are updating
      double speciesCharge;
  };
}

#endif //  LC_POISSON_BRACKET_GYRO_EQUATION_4D_H
