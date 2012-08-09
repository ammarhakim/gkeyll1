/**
 * @file	LcHyperEquation.h
 *
 * @brief	Interface to hyperbolic equations.
 */
#ifndef LC_HYPER_EQUATION_H
#define LC_HYPER_EQUATION_H

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

namespace Lucee
{
/**
 * Represents a system of hyperbolic equations. This class provides a
 * rich interface to compute various quantities (fluxes, speeds,
 * waves, eigensystem, flux Jacobians, ...) for use in various
 * numercial schemes. Not all methods are required for every scheme
 * and the documentation for a particular scheme should be consulted
 * to determine which methods are actually used.
 */
  class HyperEquation : public Lucee::BasicObj
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create a new hyperbolic equation system.
 *
 * @param Number of equations in system.
 * @param mwave Number of waves in system.
 */
      HyperEquation(unsigned meqn, unsigned mwave);

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Return number of equations in system.
 *
 * @return number of equations.
 */
      unsigned getNumEqns() const { return meqn; }

/**
 * Return number of waves in system.
 *
 * @return number of waves.
 */
      unsigned getNumWaves() const { return mwave; }

/**
 * Rotate data to local coordinate system.
 *
 * @param c Coordinate system to rotate data to.
 * @param inQ Input conserved variables.
 * @param outQ Rotated conserved variables. 
 */
      virtual void rotateToLocal(const Lucee::RectCoordSys& c, const double* inQ, double* outQ) = 0;

/**
 * Rotate data to global coordinate system.
 *
 * @param c Coordinate system to rotate data to.
 * @param inQ Input conserved variables.
 * @param outQ Rotated conserved variables. 
 */
      virtual void rotateToGlobal(const Lucee::RectCoordSys& c, const double* inQ, double* outQ) = 0;

/**
 * Compute flux for this equation system.
 *
 * @param c Coordinate system in which to compute flux.
 * @param q Conserved variables for which to compute flux.
 * @param f On output, this contains the flux.
 */
      virtual void flux(const Lucee::RectCoordSys& c, const double* q, double* f);

/**
 * Compute the minimum and maximum wave speeds in the system. s[0] is
 * the minimum wave speed and s[1] is the maximum wave speed.
 *
 * @param c Coordinate system in which to compute speeds.
 * @param q Conserved variables for which to compute speeds.
 * @param s On output, s[0] is the minimum speed and s[1] the maximum speed.
 */
      virtual void speeds(const Lucee::RectCoordSys& c, const double* q, double s[2]);

/**
 * Compute primitive variables given conserved variables.
 *
 * @param q Conserved variables for which to compute primitive variables.
 * @param v On output, primitive variables.
 */
      virtual void primitive(const double* q, double* v) const;

/**
 * Compute conserved variables given primitive variables.
 *
 * @param v Primitive variables for which to compute conserved variables.
 * @param q On output, conserved variables.
 */
      virtual void conserved(const double* v, double* q) const;

/**
 * Decompose jump into waves and wave-speeds using right and left
 * states. The states and jump are already in the local coordinate
 * system specified by 'c'. Hence, in most case (equation system is
 * isotropic) the coordinate system should be ignored.
 *
 * @param c Coordinate system in which to compute waves.
 * @param jump Jump to decompose.
 * @param ql Left state conserved variables.
 * @param qr Right state conserved variables.
 * @param waves On output, waves. This matrix has shape (meqn X mwave).
 * @param s On output, wave speeds.
 */
      virtual void waves(const Lucee::RectCoordSys& c,
        const Lucee::ConstFieldPtr<double>& jump,
        const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
        Lucee::Matrix<double>& waves, Lucee::FieldPtr<double>& s);

/**
 * Compute numerical flux for this equation system. Numerical flux
 * depends on left and right states. This method should also return
 * the maximum wave speed computed from the states. The states are
 * already in the normal-tangent space and so for isotropic systems
 * the coordinate systems can be ignored.
 *
 * @param c Coordinate system in which to compute flux.
 * @param ql Left conserved variable state.
 * @param qr Right conserved variable state.
 * @param f On output, this contains the numerical flux.
 * @return Maximum wave speed from left/right state.
 */
      virtual double numericalFlux(const Lucee::RectCoordSys& c,
        const double* ql, const double* qr, double* f);

/**
 * Project given vector on left-eigenvectors of flux-jacobian. The
 * conserved variables are already in the local coordinate system and
 * hence for isotropic equations the coordinate system 'c' can be
 * ignored.
 *
 * @param c Coordinate system in which to compute eigensystem.
 * @param q Conserved variables at which flux Jacobian is to be computed.
 * @param vec Vector to project.
 * @param coeff On output, the projection of 'vec' on left-eigenvectors.
 */
      virtual void projectOnLeftEigenvectors(const Lucee::RectCoordSys& c,
        const double* q, const double* vec, double* coeff);

/**
 * Reconstruct vector by weighted sum of right eigenvectors. The
 * conserved variables are already in the local coordinate system and
 * hence for isotropic equations the coordinate system 'c' can be
 * ignored.
 *
 * @param c Coordinate system in which to compute eigensystem.
 * @param q Conserved variables at which flux Jacobian is to be computed.
 * @param coeff Coefficients to multiply corresponding right-eigenvectors.
 * @param vec On output, the reconstructured vector.
 */
      virtual void reconWithRightEigenvectors(const Lucee::RectCoordSys& c,
        const double* q, const double* coeff, double* vec);

/**
 * Compute eigensystem for equations give a state. This method should
 * return all eigenvalues and right and left eigenvectors.
 *
 * @param c Coordinate system in which to compute eigensystem.
 * @param q State at which to compute eigensystem.
 * @param ev On output, eigenvalues of system.
 * @param rev On output, right eigenvectors of system stored as columns.
 * @param lev On output, left eigenvectors of system stored as rows.
 */
      virtual void eigensystem(const Lucee::RectCoordSys& c,
        const double *q,
        Lucee::Vector<double>& ev, Lucee::Matrix<double>& rev, Lucee::Matrix<double>& lev);

/**
 * Compute the quasi-Linear matrix given the primitive state.
 *
 * @param c Coordinate system in which to compute eigensystem.
 * @param v Primitive variables (state) at which quasi-Linear matrix.
 * @param qlMat Quasi-Linear matrix computed at 'v'.
 */
      virtual void quasiLinearMatrix(const Lucee::RectCoordSys& c,
        const double *v, Lucee::Matrix<double>& qlMat);

/**
 * Check if conserved variables satisfies invariant domains of the
 * system. Return true if it does, false otherwise.
 *
 * @param q Conserved variables.
 * @return true if invariant domains are satisfied, false otherwise.
 */
      virtual bool isInvariantDomain(const double* q) const;

/**
 * Compute fluctuations using q-waves from waves and speeds.
 *
 * @param waves Waves. This matrix has shape (meqn X mwave).
 * @param s Wave speeds.
 * @param amdq On output, fluctuations in the negative direction.
 * @param apdq On output, fluctuations in the positive direction.
 */
      void qFluctuations(const Lucee::Matrix<double>& waves, const Lucee::FieldPtr<double>& s,
        Lucee::FieldPtr<double>& amdq, Lucee::FieldPtr<double>& apdq);

/**
 * Compute fluctuations using f-waves from waves and speeds.
 *
 * @param waves Waves. This matrix has shape (meqn X mwave).
 * @param s Wave speeds.
 * @param amdq On output, fluctuations in the negative direction.
 * @param apdq On output, fluctuations in the positive direction.
 */
      void fFluctuations(const Lucee::Matrix<double>& waves, const Lucee::FieldPtr<double>& s,
        Lucee::FieldPtr<double>& amdq, Lucee::FieldPtr<double>& apdq);

/**
 * Method that performs registration of Lua functions.
 *
 * @param lfm Lua function map object.
 */
      static void appendLuaCallableMethods(Lucee::LuaFuncMap& lfm);

    protected:

    private:
/** Number of equations */
      unsigned meqn;
/** Number of waves */
      unsigned mwave;

/**
 * Lua callable method to compute primitive variables from conserved
 * variables.
 *
 * @param L Lua state to use.
 * @return number of output parameters.
 */
      static int luaPrimitive(lua_State *L);

/**
 * Lua callable method to compute conserved variables from primitive
 * variables.
 *
 * @param L Lua state to use.
 * @return number of output parameters.
 */
      static int luaConserved(lua_State *L);

/**
 * Compute primitive variables given conserved variables.
 *
 * @param cons Conserved variables to use.
 * @param prim On ouput this contains primitive variables.
 */
      template <unsigned NDIM>
      void
      calcPrimVars(const Lucee::StructGridField<NDIM, double>& cons, Lucee::StructGridField<NDIM, double>& prim) const 
      {
        Lucee::FieldPtr<double> pPtr = prim.createPtr();
        Lucee::ConstFieldPtr<double> cPtr = cons.createConstPtr();

        int idx[NDIM];
// loop over extended region and compute primitive variables
        Lucee::RowMajorSequencer<NDIM> seq(prim.getExtRegion());
        while (seq.step())
        {
          seq.fillWithIndex(idx);
          cons.setPtr(cPtr, idx);
          prim.setPtr(pPtr, idx);
// compute primitive variables
          this->primitive(&cPtr[0], &pPtr[0]);
        }
      }

/**
 * Compute conserved variables given primitive variables.
 *
 * @param prim Primitive variables to use.
 * @param cons On ouput this contains conserved variables.
 */
      template <unsigned NDIM>
      void
      calcConsVars(const Lucee::StructGridField<NDIM, double>& prim, Lucee::StructGridField<NDIM, double>& cons) const 
      {
        Lucee::ConstFieldPtr<double> pPtr = prim.createConstPtr();
        Lucee::FieldPtr<double> cPtr = cons.createPtr();

        int idx[NDIM];
// loop over extended region and compute primitive variables
        Lucee::RowMajorSequencer<NDIM> seq(cons.getExtRegion());
        while (seq.step())
        {
          seq.fillWithIndex(idx);
          cons.setPtr(cPtr, idx);
          prim.setPtr(pPtr, idx);
// compute primitive variables
          this->conserved(&pPtr[0], &cPtr[0]);
        }
      }
  };
}

#endif //  LC_HYPER_EQUATION_H
