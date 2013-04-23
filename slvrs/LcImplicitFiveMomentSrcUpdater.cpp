/**
 * @file	LcImplicitFiveMomentSrcUpdater.cpp
 *
 * @brief	Implicit updater for 5-moment source terms
 */

// lucee includes
#include <LcImplicitFiveMomentSrcUpdater.h>
#include <LcStructGridField.h>

// eigen inlcudes
#include <Eigen/Eigen>

namespace Lucee
{

  static const unsigned X = 0;
  static const unsigned Y = 1;
  static const unsigned Z = 2;

  static const unsigned RHO = 0;
  static const unsigned RHOUX = 1;
  static const unsigned RHOUY = 2;
  static const unsigned RHOUZ = 3;
  static const unsigned ER = 4;

  static const unsigned EX = 0;
  static const unsigned EY = 1;
  static const unsigned EZ = 2;
  static const unsigned BX = 3;
  static const unsigned BY = 4;
  static const unsigned BZ = 5;

  static const int COL_PIV_HOUSEHOLDER_QR = 0;
  static const int PARTIAL_PIV_LU = 1;

// set ids for module system
  template <> const char *ImplicitFiveMomentSrcUpdater<1>::id = "ImplicitFiveMomentSrc1D";
  template <> const char *ImplicitFiveMomentSrcUpdater<2>::id = "ImplicitFiveMomentSrc2D";
  template <> const char *ImplicitFiveMomentSrcUpdater<3>::id = "ImplicitFiveMomentSrc3D";

  template <unsigned NDIM>
  void
  ImplicitFiveMomentSrcUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    nFluids = (int) tbl.getNumber("numFluids");

// charge on each fluid species
    charge = tbl.getNumVec("charge");
    if (charge.size() != nFluids)
    {
      Lucee::Except lce("ImplicitFiveMomentSrcUpdater::readInput: 'charge' table should have exactly ");
      lce << nFluids << " entries. Provided " << charge.size() << " instead.";
      throw lce;
    }

// mass of each fluid species
    mass = tbl.getNumVec("mass");
    if (mass.size() != nFluids)
    {
      Lucee::Except lce("ImplicitFiveMomentSrcUpdater::readInput: 'mass' table should have exactly ");
      lce << nFluids << " entries. Provided " << mass.size() << " instead.";
      throw lce;
    }

// permittivity of free space
    epsilon0 = tbl.getNumber("epsilon0");

// linear solver type
    linSolType = PARTIAL_PIV_LU;
    if (tbl.hasString("linearSolver"))
    {
      if (tbl.getString("linearSolver") == "partialPivLu")
        linSolType = PARTIAL_PIV_LU;
      else if (tbl.getString("linearSolver") == "colPivHouseholderQr")
        linSolType = COL_PIV_HOUSEHOLDER_QR;
      else
      {
        Lucee::Except lce("ImplicitFiveMomentSrcUpdater::readInput: 'linearSolver' must be ");
        lce << " one of partialPivLu or colPivHouseholderQr. Provided " << tbl.getString("linearSolver")
            << " instead";
        throw lce;
      }
    }

    qbym.resize(nFluids);
    qbym2.resize(nFluids);
    for (unsigned i=0; i<nFluids; ++i)
    {
      qbym[i] = charge[i]/mass[i];
      qbym2[i] = qbym[i]*qbym[i];
    }
  }

  template <unsigned NDIM>
  void
  ImplicitFiveMomentSrcUpdater<NDIM>::initialize()
  {
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  ImplicitFiveMomentSrcUpdater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    double dt = t-this->getCurrTime();
    double dt1 = 0.5*dt;
    double dt2 = 0.5*dt/epsilon0;

// get pointers to each fluid
    std::vector<Lucee::Field<NDIM, double>* > fluids;
    for (unsigned i=0; i<nFluids; ++i)
      fluids.push_back(&this->getOut<Lucee::Field<NDIM, double> >(i));

// get EM field (this is last out parameter)
    Lucee::Field<NDIM, double>& emField = this->getOut<Lucee::Field<NDIM, double> >(nFluids);

    Lucee::FieldPtr<double> fPtr = fluids[0]->createPtr();
    Lucee::FieldPtr<double> emPtr = emField.createPtr();

    Eigen::MatrixXd lhs;
    lhs = Eigen::MatrixXd::Constant(3*nFluids+3, 3*nFluids+3, 0.0);
    Eigen::VectorXd rhs(3*nFluids+3);

    int idx[NDIM];
    Lucee::Region<NDIM, int> localRgn = emField.getRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);

      emField.setPtr(emPtr, idx);
// fill elements corresponding to each fluid
      for (unsigned n=0; n<nFluids; ++n)
      {
        fluids[n]->setPtr(fPtr, idx);

// eqn. for X-component of current
        lhs(fidx(n,X), fidx(n,X)) = 1.0;
        lhs(fidx(n,X), fidx(n,Y)) = -dt1*qbym[n]*emPtr[BZ];
        lhs(fidx(n,X), fidx(n,Z)) = dt1*qbym[n]*emPtr[BY];
        lhs(fidx(n,X), eidx(X)) = -dt1*qbym2[n]*fPtr[RHO];

// eqn. for Y-component of current
        lhs(fidx(n,Y), fidx(n,X)) = dt1*qbym[n]*emPtr[BZ];
        lhs(fidx(n,Y), fidx(n,Y)) = 1.0;
        lhs(fidx(n,Y), fidx(n,Z)) = -dt1*qbym[n]*emPtr[BX];
        lhs(fidx(n,Y), eidx(Y)) = -dt1*qbym2[n]*fPtr[RHO];

// eqn. for Z-component of current
        lhs(fidx(n,Z), fidx(n,X)) = -dt1*qbym[n]*emPtr[BY];
        lhs(fidx(n,Z), fidx(n,Y)) = dt1*qbym[n]*emPtr[BX];
        lhs(fidx(n,Z), fidx(n,Z)) = 1.0;
        lhs(fidx(n,Z), eidx(Z)) = -dt1*qbym2[n]*fPtr[RHO];

// fill corresponding RHS elements
        rhs(fidx(n,X)) = qbym[n]*fPtr[RHOUX];
        rhs(fidx(n,Y)) = qbym[n]*fPtr[RHOUY];
        rhs(fidx(n,Z)) = qbym[n]*fPtr[RHOUZ];

// set current contribution to electric field equation
        lhs(eidx(X), fidx(n,X)) = dt2;
        lhs(eidx(Y), fidx(n,Y)) = dt2;
        lhs(eidx(Z), fidx(n,Z)) = dt2;
      }

// fill in elements for electric field equations
      lhs(eidx(EX), eidx(EX)) = 1.0;
      lhs(eidx(EY), eidx(EY)) = 1.0;
      lhs(eidx(EZ), eidx(EZ)) = 1.0;

      rhs(eidx(EX)) = emPtr[EX];
      rhs(eidx(EY)) = emPtr[EY];
      rhs(eidx(EZ)) = emPtr[EZ];

// invert to find solution
      Eigen::VectorXd sol;
      if (linSolType == COL_PIV_HOUSEHOLDER_QR)
        sol = lhs.colPivHouseholderQr().solve(rhs);
      else if (linSolType == PARTIAL_PIV_LU)
        sol = lhs.partialPivLu().solve(rhs);
      else
      { /* can not happen */ }

      double keold = 0.0;
// update solution for fluids (solution is at half-time step)
      for (unsigned n=0; n<nFluids; ++n)
      {
        fluids[n]->setPtr(fPtr, idx);

// compute old kinetic energy before momenta are over-written
        keold = 0.5*(fPtr[RHOUX]*fPtr[RHOUX] + fPtr[RHOUY]*fPtr[RHOUY] + fPtr[RHOUZ]*fPtr[RHOUZ])/fPtr[RHO];

// momentum equation (sol has currents, so divide out charge and
// multiply by mass to give momentum density)
        fPtr[RHOUX] = 2*sol(fidx(n,X))/qbym[n] - fPtr[RHOUX];
        fPtr[RHOUY] = 2*sol(fidx(n,Y))/qbym[n] - fPtr[RHOUY];
        fPtr[RHOUZ] = 2*sol(fidx(n,Z))/qbym[n] - fPtr[RHOUZ];

// energy equation (there is no explicit energy source, so just
// recompute new kinetic energy to update total energy)
        fPtr[ER] = fPtr[ER] - keold
          + 0.5*(fPtr[RHOUX]*fPtr[RHOUX] + fPtr[RHOUY]*fPtr[RHOUY] + fPtr[RHOUZ]*fPtr[RHOUZ])/fPtr[RHO];
      }

// update electric field
      emPtr[EX] = 2*sol(eidx(X)) - emPtr[EX];
      emPtr[EY] = 2*sol(eidx(Y)) - emPtr[EY];
      emPtr[EZ] = 2*sol(eidx(Z)) - emPtr[EZ];
    }
    
    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  ImplicitFiveMomentSrcUpdater<NDIM>::declareTypes()
  {
    this->setLastOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class ImplicitFiveMomentSrcUpdater<1>;
  template class ImplicitFiveMomentSrcUpdater<2>;
  template class ImplicitFiveMomentSrcUpdater<3>;
}

