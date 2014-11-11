/**
 * @file	LcPoissonBracketGyroEquation5D.h
 *
 * @brief	Compute Jacobian*alpha at all quadrature points
 */
#ifndef LC_POISSON_BRACKET_GYRO_EQUATION_5D_H
#define LC_POISSON_BRACKET_GYRO_EQUATION_5D_H

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
  class PoissonBracketGyroEquation5D : public Lucee::PoissonBracketEquation
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create a new equation.
 */
      PoissonBracketGyroEquation5D();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

      virtual void computeAlphaAtQuadNodes(const Eigen::MatrixXd& hamiltonian, const Eigen::MatrixXd& interpMat,
          const int idx[], Eigen::MatrixXd& alpha);
    private:
      Lucee::Field<5, double> *bStarYField;
      Lucee::Field<5, double> *bStarZField;
      // Mass of species we are updating
      double speciesMass;
      // Charge of species we are updating
      double speciesCharge;
  };
}

#endif //  LC_POISSON_BRACKET_GYRO_EQUATION_5D_H
