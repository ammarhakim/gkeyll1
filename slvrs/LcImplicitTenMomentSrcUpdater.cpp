/**
 * @file	LcImplicitTenMomentSrcUpdater.cpp
 *
 * @brief	Implicit updater for 5-moment source terms
 */

// lucee includes
#include <LcImplicitTenMomentSrcUpdater.h>
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

  static const unsigned P11 = 4;
  static const unsigned P12 = 5;
  static const unsigned P13 = 6;
  static const unsigned P22 = 7;
  static const unsigned P23 = 8;
  static const unsigned P33 = 9;

  static const unsigned EX = 0;
  static const unsigned EY = 1;
  static const unsigned EZ = 2;
  static const unsigned BX = 3;
  static const unsigned BY = 4;
  static const unsigned BZ = 5;

  static const int COL_PIV_HOUSEHOLDER_QR = 0;
  static const int PARTIAL_PIV_LU = 1;

// set ids for module system
  template <> const char *ImplicitTenMomentSrcUpdater<1>::id = "ImplicitTenMomentSrc1D";
  template <> const char *ImplicitTenMomentSrcUpdater<2>::id = "ImplicitTenMomentSrc2D";
  template <> const char *ImplicitTenMomentSrcUpdater<3>::id = "ImplicitTenMomentSrc3D";

  template <unsigned NDIM>
  void
  ImplicitTenMomentSrcUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    nFluids = (int) tbl.getNumber("numFluids");

// charge on each fluid species
    charge = tbl.getNumVec("charge");
    if (charge.size() != nFluids)
    {
      Lucee::Except lce("ImplicitTenMomentSrcUpdater::readInput: 'charge' table should have exactly ");
      lce << nFluids << " entries. Provided " << charge.size() << " instead.";
      throw lce;
    }

// mass of each fluid species
    mass = tbl.getNumVec("mass");
    if (mass.size() != nFluids)
    {
      Lucee::Except lce("ImplicitTenMomentSrcUpdater::readInput: 'mass' table should have exactly ");
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
        Lucee::Except lce("ImplicitTenMomentSrcUpdater::readInput: 'linearSolver' must be ");
        lce << " one of partialPivLu or colPivHouseholderQr. Provided " << tbl.getString("linearSolver")
            << " instead";
        throw lce;
      }
    }

// flag to indicate if there is a static field
    hasStatic = false;
    if (tbl.hasBool("hasStaticField"))
      hasStatic = tbl.getBool("hasStaticField");

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
  ImplicitTenMomentSrcUpdater<NDIM>::initialize()
  {
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  ImplicitTenMomentSrcUpdater<NDIM>::update(double t)
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

// get static EM field if one is present (it is the first and only
// input field)
    const Lucee::Field<NDIM, double>* staticField = 0;
    if (hasStatic)
      staticField = &this->getInp<Lucee::Field<NDIM, double> >(0);

    Lucee::FieldPtr<double> fPtr = fluids[0]->createPtr();
    Lucee::FieldPtr<double> emPtr = emField.createPtr();

    std::vector<double> zeros(6);
    for (unsigned i=0; i<6; ++i) zeros[i] = 0.0;
    Lucee::ConstFieldPtr<double> staticEmPtr(zeros);

// for momentum and electric field update
    Eigen::MatrixXd lhs;
    lhs = Eigen::MatrixXd::Constant(3*nFluids+3, 3*nFluids+3, 0.0);
    Eigen::VectorXd rhs(3*nFluids+3);

// for pressure update
    Eigen::MatrixXd prLhs;
    prLhs = Eigen::MatrixXd::Constant(6, 6, 0.0);
    Eigen::VectorXd prRhs(6);
// updated pressure tensor
    std::vector<double> prTen(6*nFluids);

    int idx[NDIM];
    Lucee::Region<NDIM, int> localRgn = emField.getRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);

      emField.setPtr(emPtr, idx);
      if (hasStatic)
        staticField->setPtr(staticEmPtr, idx);

// fill elements corresponding to each fluid
      for (unsigned n=0; n<nFluids; ++n)
      {
        fluids[n]->setPtr(fPtr, idx);

// eqn. for X-component of current
        lhs(fidx(n,X), fidx(n,X)) = 1.0;
        lhs(fidx(n,X), fidx(n,Y)) = -dt1*qbym[n]*(emPtr[BZ]+staticEmPtr[BZ]);
        lhs(fidx(n,X), fidx(n,Z)) = dt1*qbym[n]*(emPtr[BY]+staticEmPtr[BY]);
        lhs(fidx(n,X), eidx(X)) = -dt1*qbym2[n]*fPtr[RHO];

// eqn. for Y-component of current
        lhs(fidx(n,Y), fidx(n,X)) = dt1*qbym[n]*(emPtr[BZ]+staticEmPtr[BZ]);
        lhs(fidx(n,Y), fidx(n,Y)) = 1.0;
        lhs(fidx(n,Y), fidx(n,Z)) = -dt1*qbym[n]*(emPtr[BX]+staticEmPtr[BX]);
        lhs(fidx(n,Y), eidx(Y)) = -dt1*qbym2[n]*fPtr[RHO];

// eqn. for Z-component of current
        lhs(fidx(n,Z), fidx(n,X)) = -dt1*qbym[n]*(emPtr[BY]+staticEmPtr[BY]);
        lhs(fidx(n,Z), fidx(n,Y)) = dt1*qbym[n]*(emPtr[BX]+staticEmPtr[BX]);
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

// compute increments from pressure tensor terms
      for (unsigned n=0; n<nFluids; ++n)
      {
        fluids[n]->setPtr(fPtr, idx);

// assemble LHS terms
        prLhs(0,0) = 1;
        prLhs(0,1) = -2*dt1*qbym[n]*emPtr[BZ];
        prLhs(0,2) = 2*dt1*qbym[n]*emPtr[BY];
        prLhs(0,3) = 0;
        prLhs(0,4) = 0;
        prLhs(0,5) = 0;
        prLhs(1,0) = dt1*qbym[n]*emPtr[BZ];
        prLhs(1,1) = 1;
        prLhs(1,2) = -dt1*qbym[n]*emPtr[BX];
        prLhs(1,3) = -dt1*qbym[n]*emPtr[BZ];
        prLhs(1,4) = dt1*qbym[n]*emPtr[BY];
        prLhs(1,5) = 0;
        prLhs(2,0) = -dt1*qbym[n]*emPtr[BY];
        prLhs(2,1) = dt1*qbym[n]*emPtr[BX];
        prLhs(2,2) = 1;
        prLhs(2,3) = 0;
        prLhs(2,4) = -dt1*qbym[n]*emPtr[BZ];
        prLhs(2,5) = dt1*qbym[n]*emPtr[BY];
        prLhs(3,0) = 0;
        prLhs(3,1) = 2*dt1*qbym[n]*emPtr[BZ];
        prLhs(3,2) = 0;
        prLhs(3,3) = 1;
        prLhs(3,4) = -2*dt1*qbym[n]*emPtr[BX];
        prLhs(3,5) = 0;
        prLhs(4,0) = 0;
        prLhs(4,1) = -dt1*qbym[n]*emPtr[BY];
        prLhs(4,2) = dt1*qbym[n]*emPtr[BZ];
        prLhs(4,3) = dt1*qbym[n]*emPtr[BX];
        prLhs(4,4) = 1;
        prLhs(4,5) = -dt1*qbym[n]*emPtr[BX];
        prLhs(5,0) = 0;
        prLhs(5,1) = 0;
        prLhs(5,2) = -2*dt1*qbym[n]*emPtr[BY];
        prLhs(5,3) = 0;
        prLhs(5,4) = 2*dt1*qbym[n]*emPtr[BX];
        prLhs(5,5) = 1;

// RHS matrix
        prRhs[0] = fPtr[P11] - fPtr[RHOUX]*fPtr[RHOUX]/fPtr[RHO];
        prRhs[1] = fPtr[P12] - fPtr[RHOUX]*fPtr[RHOUY]/fPtr[RHO];
        prRhs[2] = fPtr[P13] - fPtr[RHOUX]*fPtr[RHOUZ]/fPtr[RHO];
        prRhs[3] = fPtr[P22] - fPtr[RHOUY]*fPtr[RHOUY]/fPtr[RHO];
        prRhs[4] = fPtr[P23] - fPtr[RHOUY]*fPtr[RHOUZ]/fPtr[RHO];
        prRhs[5] = fPtr[P33] - fPtr[RHOUZ]*fPtr[RHOUZ]/fPtr[RHO];

// solve to compute increment in pressure tensor
        Eigen::VectorXd prSol;
        if (linSolType == COL_PIV_HOUSEHOLDER_QR)
          prSol = prLhs.colPivHouseholderQr().solve(prRhs);
        else if (linSolType == PARTIAL_PIV_LU)
          prSol = prLhs.partialPivLu().solve(prRhs);
        else
        { /* can not happen */ }

// update solution
        for (unsigned i=0; i<6; ++i)
          prTen[6*n+i] = 2*prSol[i]-prRhs[i];
      }

// update solution for fluids
      for (unsigned n=0; n<nFluids; ++n)
      {
        fluids[n]->setPtr(fPtr, idx);

// momentum equation (sol has currents, so divide out charge and
// multiply by mass to give momentum density)
        fPtr[RHOUX] = 2*sol(fidx(n,X))/qbym[n] - fPtr[RHOUX];
        fPtr[RHOUY] = 2*sol(fidx(n,Y))/qbym[n] - fPtr[RHOUY];
        fPtr[RHOUZ] = 2*sol(fidx(n,Z))/qbym[n] - fPtr[RHOUZ];

// total pressure tensor equations (we have computed pressure tensor,
// so we need to add in the Reynold stress terms to get the total
// pressure)
        fPtr[P11] = fPtr[RHOUX]*fPtr[RHOUX]/fPtr[RHO] + prTen[6*n+0];
        fPtr[P12] = fPtr[RHOUX]*fPtr[RHOUY]/fPtr[RHO] + prTen[6*n+1];
        fPtr[P13] = fPtr[RHOUX]*fPtr[RHOUZ]/fPtr[RHO] + prTen[6*n+2];
        fPtr[P22] = fPtr[RHOUY]*fPtr[RHOUY]/fPtr[RHO] + prTen[6*n+3];
        fPtr[P23] = fPtr[RHOUY]*fPtr[RHOUZ]/fPtr[RHO] + prTen[6*n+4];
        fPtr[P33] = fPtr[RHOUZ]*fPtr[RHOUZ]/fPtr[RHO] + prTen[6*n+5];
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
  ImplicitTenMomentSrcUpdater<NDIM>::declareTypes()
  {
    this->setLastOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class ImplicitTenMomentSrcUpdater<1>;
  template class ImplicitTenMomentSrcUpdater<2>;
  template class ImplicitTenMomentSrcUpdater<3>;
}

