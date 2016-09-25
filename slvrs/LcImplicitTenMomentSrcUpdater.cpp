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
  static const int ANALYTIC = 2;

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
      else if (tbl.getString("linearSolver") == "analytic")
        linSolType = ANALYTIC;
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

      // variables needed for analytic solver
      std::vector<double> lambda(nFluids);
      std::vector<double> omega(nFluids);
      double gamma2 = 0.0;
      double delta = 0.0;
      double theta = 0.0; 
      double bx = (emPtr[BX] + staticEmPtr[BX]);
      double by = (emPtr[BY] + staticEmPtr[BY]);
      double bz = (emPtr[BZ] + staticEmPtr[BZ]);
      Eigen::Vector3d B; 
      Eigen::Vector3d E;
      Eigen::Vector3d F(0,0,0); 
      Eigen::Vector3d K(0,0,0);// the k vector used to update the implicit solution in Smithe(2007)
      B(0) = bx; B(1) = by; B(2) = bz;
      E(0) = emPtr[EX]; E(1) = emPtr[EY]; E(2) = emPtr[EZ];
      double babs = std::sqrt(bx*bx + by*by + bz*bz);


// fill elements corresponding to each fluid
      for (unsigned n=0; n<nFluids; ++n)
      {
        fluids[n]->setPtr(fPtr, idx);
        
        if (linSolType != ANALYTIC) {
// eqn. for X-component of current
          lhs(fidx(n,X), fidx(n,X)) = 1.0;
          lhs(fidx(n,X), fidx(n,Y)) = -dt1*qbym[n]*bz;
          lhs(fidx(n,X), fidx(n,Z)) = dt1*qbym[n]*by;
          lhs(fidx(n,X), eidx(X)) = -dt1*qbym2[n]*fPtr[RHO];

// eqn. for Y-component of current
          lhs(fidx(n,Y), fidx(n,X)) = dt1*qbym[n]*bz;
          lhs(fidx(n,Y), fidx(n,Y)) = 1.0;
          lhs(fidx(n,Y), fidx(n,Z)) = -dt1*qbym[n]*bx;
          lhs(fidx(n,Y), eidx(Y)) = -dt1*qbym2[n]*fPtr[RHO];

// eqn. for Z-component of current
          lhs(fidx(n,Z), fidx(n,X)) = -dt1*qbym[n]*by;
          lhs(fidx(n,Z), fidx(n,Y)) = dt1*qbym[n]*bx;
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
        } else {
          lambda[n] = 1.00; //collisional parameter
          omega[n] = qbym[n]*dt; // omega/B
          double wp2 = fPtr[RHO]*qbym2[n]/epsilon0*dt*dt;
          //          wp2 += fPtr[RHO]*qbym2[n]/epsilon0*dt*dt;
          gamma2 += wp2*qbym2[n]*dt*dt/(1+0.25*qbym2[n]*babs*babs*dt*dt);
          delta += wp2*qbym[n]*dt/(1+0.25*qbym2[n]*babs*babs*dt*dt);
          theta += wp2/(1+qbym2[n]*dt*dt/4*babs*babs);
          Eigen::Vector3d j(fPtr[RHOUX]*qbym[n], fPtr[RHOUY]*qbym[n], fPtr[RHOUZ]*qbym[n]);

          K -= dt/2.0 *(1.0 + lambda[n])*(1.0*j + 0.25*omega[n]*omega[n]*B * B.dot(j) - 0.5*omega[n]*B.cross(j) )/(1.0+0.25*omega[n]*omega[n]*babs*babs);
          F -= dt/2.0*(1/4.0*omega[n]*omega[n]*fPtr[RHO]*qbym2[n]/epsilon0*dt/2.0*B*B.dot(E) + 1/2.0 * fPtr[RHO]*qbym2[n]/epsilon0*dt*E - 1/2.0*omega[n]*fPtr[RHO]*qbym2[n]/epsilon0*dt*B.cross(E)/2.0 )/(1.0+0.25*omega[n]*omega[n]*babs*babs);
        }
      }
      F = F + E;
// fill in elements for electric field equations
      Eigen::VectorXd sol;
      if (linSolType != ANALYTIC){
        lhs(eidx(EX), eidx(EX)) = 1.0;
        lhs(eidx(EY), eidx(EY)) = 1.0;
        lhs(eidx(EZ), eidx(EZ)) = 1.0;
        
        rhs(eidx(EX)) = emPtr[EX];
        rhs(eidx(EY)) = emPtr[EY];
        rhs(eidx(EZ)) = emPtr[EZ];

// invert to find solution
        if (linSolType == COL_PIV_HOUSEHOLDER_QR)
          sol = lhs.colPivHouseholderQr().solve(rhs);
        else if (linSolType == PARTIAL_PIV_LU)
          sol = lhs.partialPivLu().solve(rhs);
        else
          { /* can not happen */ }
      } else {

        double denom = ((1+theta/4)*(1+theta/4) + delta*delta*babs*babs/64.0) * (1 + theta/4 + gamma2/16*babs*babs);        
        double coeff_e = ((1+theta/4)*(1+theta/4) + gamma2*babs*babs/16*(1+theta/4));
        double coeff_dot = (1.0/64.0*delta*delta - 1.0/16.0*gamma2*(1+theta/4));
        double coeff_cross = 0.125*delta*(1 + gamma2*babs*babs/16 + theta/4);
        Eigen::Vector3d fk = 1.0*F + K/epsilon0; // vector holding sum of K + 2 epsilon0 E used in many calculations
        emPtr[EX] = (coeff_e*fk[EX] + coeff_dot*B(0)*B.dot(fk) + coeff_cross*(B.cross(fk)(0)) )/denom; 
        emPtr[EY] = (coeff_e*fk[EY] + coeff_dot*B(1)*B.dot(fk) + coeff_cross*(B.cross(fk)(1)) )/denom; 
        emPtr[EZ] = (coeff_e*fk[EZ] + coeff_dot*B(2)*B.dot(fk) + coeff_cross*(B.cross(fk)(2)) )/denom; 
        // update the stored E field to E_n+1/2
        E(0) = (emPtr[EX] + E(0))/2.0;
        E(1) = (emPtr[EY] + E(1))/2.0;
        E(2) = (emPtr[EZ] + E(2))/2.0;

      }
// compute increments from pressure tensor terms
      for (unsigned n=0; n<nFluids; ++n)
      {
        fluids[n]->setPtr(fPtr, idx);

        if (linSolType != ANALYTIC) {
// assemble LHS terms
        prLhs(0,0) = 1;
        prLhs(0,1) = -2*dt1*qbym[n]*bz;
        prLhs(0,2) = 2*dt1*qbym[n]*by;
        prLhs(0,3) = 0;
        prLhs(0,4) = 0;
        prLhs(0,5) = 0;
        prLhs(1,0) = dt1*qbym[n]*bz;
        prLhs(1,1) = 1;
        prLhs(1,2) = -dt1*qbym[n]*bx;
        prLhs(1,3) = -dt1*qbym[n]*bz;
        prLhs(1,4) = dt1*qbym[n]*by;
        prLhs(1,5) = 0;
        prLhs(2,0) = -dt1*qbym[n]*by;
        prLhs(2,1) = dt1*qbym[n]*bx;
        prLhs(2,2) = 1;
        prLhs(2,3) = 0;
        prLhs(2,4) = -dt1*qbym[n]*bz;
        prLhs(2,5) = dt1*qbym[n]*by;
        prLhs(3,0) = 0;
        prLhs(3,1) = 2*dt1*qbym[n]*bz;
        prLhs(3,2) = 0;
        prLhs(3,3) = 1;
        prLhs(3,4) = -2*dt1*qbym[n]*bx;
        prLhs(3,5) = 0;
        prLhs(4,0) = 0;
        prLhs(4,1) = -dt1*qbym[n]*by;
        prLhs(4,2) = dt1*qbym[n]*bz;
        prLhs(4,3) = dt1*qbym[n]*bx;
        prLhs(4,4) = 1;
        prLhs(4,5) = -dt1*qbym[n]*bx;
        prLhs(5,0) = 0;
        prLhs(5,1) = 0;
        prLhs(5,2) = -2*dt1*qbym[n]*by;
        prLhs(5,3) = 0;
        prLhs(5,4) = 2*dt1*qbym[n]*bx;
        prLhs(5,5) = 1;
        }
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
        else if (linSolType == ANALYTIC) {
          // useful variables
          double dtsq = dt1*dt1;
          double dt3 = dtsq*dt1;
          double dt4 = dt3*dt1;
          double qb2 = qbym[n]*qbym[n];
          double qb3 = qbym[n]*qb2;
          double qb4 = qbym[n]*qb3;
          double Bx2 = bx*bx;
          double Bx3 = bx*Bx2;
          double Bx4 = bx*Bx3;
          double By2 = by*by;
          double By3 = by*By2;
          double By4 = by*By3;
          double Bz2 = bz*bz;
          double Bz3 = bz*Bz2;
          double Bz4 = bz*Bz3;
          double d = 1 + 5*(Bx2 + By2 + Bz2)*dtsq*qb2 + 4*(Bx2 + By2 + Bz2)*(Bx2 + By2 + Bz2)*dt4*qb4;
          prSol = Eigen::VectorXd(6);
          prSol[0] = (prRhs[0] + 2*dt1*(bz*prRhs[1] - by*prRhs[2])*qbym[n] + dtsq*(5*Bx2*prRhs[0] + 2*bx*(by*prRhs[1] + bz*prRhs[2]) + Bz2*(3*prRhs[0] + 2*prRhs[3]) - 4*by*bz*prRhs[4] + By2*(3*prRhs[0] + 2*prRhs[5]))*qb2 + 2*dt3*(4*Bx2*(bz*prRhs[1] - by*prRhs[2]) - (By2 + Bz2)*(-(bz*prRhs[1]) + by*prRhs[2]) - 3*bx*(By2*prRhs[4] - Bz2*prRhs[4] + by*bz*(-prRhs[3] + prRhs[5])))*qb3 + 2*dt4*(2*Bx4*prRhs[0] + 4*Bx3*(by*prRhs[1] + bz*prRhs[2]) - 2*bx*(By2 + Bz2)*(by*prRhs[1] + bz*prRhs[2]) + (By2 + Bz2)*(Bz2*(prRhs[0] + prRhs[3]) - 2*by*bz*prRhs[4] + By2*(prRhs[0] + prRhs[5])) + Bx2*(4*by*bz*prRhs[4] + By2*(3*prRhs[3] + prRhs[5]) + Bz2*(prRhs[3] + 3*prRhs[5])))*qb4)/d;
          prSol[1] =  (prRhs[1] + dt1*(bx*prRhs[2] + bz*(-prRhs[0] + prRhs[3]) - by*prRhs[4])*qbym[n] + dtsq*(4*Bx2*prRhs[1] + 4*By2*prRhs[1] + Bz2*prRhs[1] + 3*by*bz*prRhs[2] + bx*(3*bz*prRhs[4] + by*(prRhs[0] + prRhs[3] - 2*prRhs[5])))*qb2 + dt3*(4*Bx3*prRhs[2] - 2*bx*(By2 + Bz2)*prRhs[2] + Bz3*(-prRhs[0] + prRhs[3]) - 4*By3*prRhs[4] + 2*by*Bz2*prRhs[4] - By2*bz*(prRhs[0] - 4*prRhs[3] + 3*prRhs[5]) + Bx2*(2*by*prRhs[4] + bz*(-4*prRhs[0] + prRhs[3] + 3*prRhs[5])))*qb3 + 2*bx*by*dt4*(6*bx*(by*prRhs[1] + bz*prRhs[2]) + 6*by*bz*prRhs[4] - Bz2*(prRhs[0] + prRhs[3] - 2*prRhs[5]) + Bx2*(2*prRhs[0] - prRhs[3] - prRhs[5]) - By2*(prRhs[0] - 2*prRhs[3] + prRhs[5]))*qb4)/d;
          prSol[2] =  (prRhs[2] + dt1*(-(bx*prRhs[1]) + bz*prRhs[4] + by*(prRhs[0] - prRhs[5]))*qbym[n] + dtsq*(3*by*bz*prRhs[1] + 4*Bx2*prRhs[2] + By2*prRhs[2] + 4*Bz2*prRhs[2] + bx*(3*by*prRhs[4] + bz*(prRhs[0] - 2*prRhs[3] + prRhs[5])))*qb2 + dt3*(-4*Bx3*prRhs[1] + 2*bx*(By2 + Bz2)*prRhs[1] - 2*By2*bz*prRhs[4] + 4*Bz3*prRhs[4] + by*Bz2*(prRhs[0] + 3*prRhs[3] - 4*prRhs[5]) + By3*(prRhs[0] - prRhs[5]) - Bx2*(2*bz*prRhs[4] + by*(-4*prRhs[0] + 3*prRhs[3] + prRhs[5])))*qb3 + 2*bx*bz*dt4*(6*bx*(by*prRhs[1] + bz*prRhs[2]) + 6*by*bz*prRhs[4] - Bz2*(prRhs[0] + prRhs[3] - 2*prRhs[5]) + Bx2*(2*prRhs[0] - prRhs[3] - prRhs[5]) - By2*(prRhs[0] - 2*prRhs[3] + prRhs[5]))*qb4)/d;
          prSol[3] =  (prRhs[3] + (-2*bz*dt1*prRhs[1] + 2*bx*dt1*prRhs[4])*qbym[n] + dtsq*(2*bx*by*prRhs[1] + 5*By2*prRhs[3] + Bz2*(2*prRhs[0] + 3*prRhs[3]) + bz*(-4*bx*prRhs[2] + 2*by*prRhs[4]) + Bx2*(3*prRhs[3] + 2*prRhs[5]))*qb2 + 2*dt3*(Bx2*(-(bz*prRhs[1]) + 3*by*prRhs[2]) - bz*(4*By2*prRhs[1] + Bz2*prRhs[1] + 3*by*bz*prRhs[2]) + Bx3*prRhs[4] + bx*(4*By2*prRhs[4] + Bz2*prRhs[4] + 3*by*bz*(-prRhs[0] + prRhs[5])))*qb3 + 2*dt4*(-2*Bx3*(by*prRhs[1] + bz*prRhs[2]) + 2*bx*(2*By2 - Bz2)*(by*prRhs[1] + bz*prRhs[2]) + 2*By4*prRhs[3] + Bz4*(prRhs[0] + prRhs[3]) + 4*By3*bz*prRhs[4] - 2*by*Bz3*prRhs[4] + Bx4*(prRhs[3] + prRhs[5]) + By2*Bz2*(prRhs[0] + 3*prRhs[5]) + Bx2*(-2*by*bz*prRhs[4] + By2*(3*prRhs[0] + prRhs[5]) + Bz2*(prRhs[0] + 2*prRhs[3] + prRhs[5])))*qb4)/d;
          prSol[4] =  (prRhs[4] + dt1*(by*prRhs[1] - bz*prRhs[2] + bx*(-prRhs[3] + prRhs[5]))*qbym[n] + dtsq*(3*bx*bz*prRhs[1] + Bx2*prRhs[4] + 4*By2*prRhs[4] + 4*Bz2*prRhs[4] + by*(3*bx*prRhs[2] + bz*(-2*prRhs[0] + prRhs[3] + prRhs[5])))*qb2 + dt3*(4*By3*prRhs[1] - 2*by*Bz2*prRhs[1] + 2*By2*bz*prRhs[2] - 4*Bz3*prRhs[2] + Bx2*(-2*by*prRhs[1] + 2*bz*prRhs[2]) + Bx3*(-prRhs[3] + prRhs[5]) + bx*(-(Bz2*(3*prRhs[0] + prRhs[3] - 4*prRhs[5])) + By2*(3*prRhs[0] - 4*prRhs[3] + prRhs[5])))*qb3 - 2*by*bz*dt4*(-6*bx*(by*prRhs[1] + bz*prRhs[2]) - 6*by*bz*prRhs[4] + Bz2*(prRhs[0] + prRhs[3] - 2*prRhs[5]) + By2*(prRhs[0] - 2*prRhs[3] + prRhs[5]) + Bx2*(-2*prRhs[0] + prRhs[3] + prRhs[5]))*qb4)/d;
          prSol[5] =  (prRhs[5] + 2*dt1*(by*prRhs[2] - bx*prRhs[4])*qbym[n] + dtsq*(2*bx*bz*prRhs[2] + by*(-4*bx*prRhs[1] + 2*bz*prRhs[4]) + 5*Bz2*prRhs[5] + By2*(2*prRhs[0] + 3*prRhs[5]) + Bx2*(2*prRhs[3] + 3*prRhs[5]))*qb2 - 2*dt3*(Bx2*(3*bz*prRhs[1] - by*prRhs[2]) - by*(3*by*bz*prRhs[1] + By2*prRhs[2] + 4*Bz2*prRhs[2]) + Bx3*prRhs[4] + bx*(3*by*bz*(-prRhs[0] + prRhs[3]) + By2*prRhs[4] + 4*Bz2*prRhs[4]))*qb3 + 2*dt4*(-2*Bx3*(by*prRhs[1] + bz*prRhs[2]) - 2*bx*(By2 - 2*Bz2)*(by*prRhs[1] + bz*prRhs[2]) + By2*Bz2*(prRhs[0] + 3*prRhs[3]) - 2*By3*bz*prRhs[4] + 4*by*Bz3*prRhs[4] + 2*Bz4*prRhs[5] + By4*(prRhs[0] + prRhs[5]) + Bx4*(prRhs[3] + prRhs[5]) + Bx2*(Bz2*(3*prRhs[0] + prRhs[3]) - 2*by*bz*prRhs[4] + By2*(prRhs[0] + prRhs[3] + 2*prRhs[5])))*qb4)/d;
        }
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
        if (linSolType != ANALYTIC) {
          fPtr[RHOUX] = 2*sol(fidx(n,X))/qbym[n] - fPtr[RHOUX];
          fPtr[RHOUY] = 2*sol(fidx(n,Y))/qbym[n] - fPtr[RHOUY];
          fPtr[RHOUZ] = 2*sol(fidx(n,Z))/qbym[n] - fPtr[RHOUZ];
        } else {
          double j_coeff = (lambda[n] - 0.25*omega[n]*omega[n]*babs*babs)/(1.0+0.25*omega[n]*omega[n]*babs*babs);
          double f_coeff = fPtr[RHO]*qbym2[n]/epsilon0*dt/(1.0+0.25*omega[n]*omega[n]*babs*babs);
          
          double dot_coeff = 0.25*omega[n]*omega[n]/(1.0+0.25*omega[n]*omega[n]*babs*babs);
          double cross_coeff = omega[n]/2.0/(1.0+0.25*omega[n]*omega[n]*babs*babs);
          // update E with the new solution
          Eigen::Vector3d j(fPtr[RHOUX]*qbym[n],fPtr[RHOUY]*qbym[n],fPtr[RHOUZ]*qbym[n]);
          Eigen::Vector3d jf = (1.0+lambda[n])*j + E*epsilon0*dt*fPtr[RHO]*qbym2[n]/epsilon0; // There is an error in smithe's paper here: F in equation (15) must be multiplied by a factor of wp^2 dt or the units are wrong.
          Eigen::Vector3d solj = j_coeff*j + f_coeff*E*epsilon0 + dot_coeff*B*B.dot(jf) - cross_coeff*B.cross(jf); 

          fPtr[RHOUX] = solj(0)/qbym[n];
          fPtr[RHOUY] = solj(1)/qbym[n];
          fPtr[RHOUZ] = solj(2)/qbym[n];

        }

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
      if (linSolType != ANALYTIC){ 
        emPtr[EX] = 2*sol(eidx(X)) - emPtr[EX];
        emPtr[EY] = 2*sol(eidx(Y)) - emPtr[EY];
        emPtr[EZ] = 2*sol(eidx(Z)) - emPtr[EZ];
      }
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

