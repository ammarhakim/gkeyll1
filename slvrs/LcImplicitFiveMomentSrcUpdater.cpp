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
  static const unsigned PHIE = 6;
  static const unsigned PHIM = 7;

  static const int COL_PIV_HOUSEHOLDER_QR = 0;
  static const int PARTIAL_PIV_LU = 1;
  static const int ANALYTIC = 2;
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
      else if (tbl.getString("linearSolver") == "analytic") {
        linSolType = ANALYTIC;
        if (tbl.hasNumber("gravity")){
        Lucee::Except lce("ImplicitFiveMomentSrcUpdater::readInput: gravity and not implemented in analytic solver "); 
        }
      }
      else
      {
        Lucee::Except lce("ImplicitFiveMomentSrcUpdater::readInput: 'linearSolver' must be ");
        lce << " one of partialPivLu or colPivHouseholderQr. Provided " << tbl.getString("linearSolver")
            << " instead";
        throw lce;
      }
    }

// flag to indicate if there is a static field
    hasStatic = false;
    if (tbl.hasBool("hasStaticField"))
      hasStatic = tbl.getBool("hasStaticField");

// flag to indicate if there is a pressure equation
    hasPressure = true;
    if (tbl.hasBool("hasPressure"))
      hasPressure = tbl.getBool("hasPressure");

// electric and magnetic error propagation speeds
    chi_e = 0.0;
    if (tbl.hasNumber("elcErrorSpeedFactor"))
      chi_e = tbl.getNumber("elcErrorSpeedFactor");

    chi_m = 0.0;
    if (tbl.hasNumber("mgnErrorSpeedFactor"))
      chi_m = tbl.getNumber("mgnErrorSpeedFactor");

// electric and magnetic damping factor
    damp_e = 0.0;
    if (tbl.hasNumber("elcErrorDampFactor"))
      damp_e = tbl.getNumber("elcErrorDampFactor");

    damp_m = 0.0;
    if (tbl.hasNumber("mgnErrorDampFactor"))
      damp_m = tbl.getNumber("mgnErrorDampFactor");

    grvDir = 0;
    gravity = 0.0;
    if (tbl.hasNumber("gravity"))
    {
      gravity = tbl.getNumber("gravity");
      if (tbl.hasNumber("dir"))
        grvDir = (unsigned) tbl.getNumber("dir");
      else
        throw Lucee::Except
          ("ImplicitFiveMomentSrcUpdater::readInput: If \"gravity\" is specified, \"dir\" must also be specified.");
    }

    qbym.resize(nFluids);
    qbym2.resize(nFluids);
    for (unsigned i=0; i<nFluids; ++i)
    {
      qbym[i] = charge[i]/mass[i];
      qbym2[i] = qbym[i]*qbym[i];
    }
    
    if (tbl.hasBool("implicitB"))
      implicitB = tbl.getBool("implicitB");
    else 
      implicitB = false;
    
    if (implicitB && linSolType != ANALYTIC) {
      throw Lucee::Except("ImplicitFiveMomentSrcUpdater::reatInput: Implicit B field update only works with analytic solver"); 
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
    
    Lucee::Field<NDIM, double>* elfNew = 0;
    if (implicitB) { 
      elfNew = &this->getOut<Lucee::Field<NDIM, double> >(nFluids+1);
    }

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

    std::vector<double> zeros3(3);
    for (unsigned i=0; i<3;++i) zeros3[i] = 0.0;
    Lucee::FieldPtr<double> elfNewPtr(zeros3);

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
      if (hasStatic)
        staticField->setPtr(staticEmPtr, idx);
      if (implicitB)
        elfNew->setPtr(elfNewPtr, idx);
      std::vector<double> lambda(nFluids);
      std::vector<double> omega(nFluids);
      //      std::vector<double> K(3); 
      //      double wp2 = 0.0;
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
      Eigen::Vector3d curlBdt(0,0,0);
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
          
// add in gravity source term to appropriate current equation
          rhs(fidx(n,grvDir)) += qbym[n]*fPtr[RHO]*gravity*dt1;
          
          // set current contribution to electric field equation
          lhs(eidx(X), fidx(n,X)) = dt2;
          lhs(eidx(Y), fidx(n,Y)) = dt2;
          lhs(eidx(Z), fidx(n,Z)) = dt2;
        } else {
          // use the analytic solver developed by Smithe (2007)
          lambda[n] = 1.00; //collisional parameter
          omega[n] = qbym[n]*dt; // omega/B
          double wp2 = fPtr[RHO]*qbym2[n]/epsilon0*dt*dt;
          //          wp2 += fPtr[RHO]*qbym2[n]/epsilon0*dt*dt;
          gamma2 += wp2*qbym2[n]*dt*dt/(1+0.25*qbym2[n]*babs*babs*dt*dt);
          delta += wp2*qbym[n]*dt/(1+0.25*qbym2[n]*babs*babs*dt*dt);
          theta += wp2/(1+qbym2[n]*dt*dt/4*babs*babs);
          Eigen::Vector3d j(fPtr[RHOUX]*qbym[n], fPtr[RHOUY]*qbym[n], fPtr[RHOUZ]*qbym[n]);
          double bdotj = B.dot(j);
          Eigen::Vector3d bcrossj = B.cross(j);

          K -= dt/2.0 *(1.0 + lambda[n])*(1.0*j + 0.25*omega[n]*omega[n]*B * bdotj - 0.5*omega[n]*bcrossj )/(1.0+0.25*omega[n]*omega[n]*babs*babs);
          F -= dt/2.0*(1/4.0*omega[n]*omega[n]*fPtr[RHO]*qbym2[n]/epsilon0*dt/2.0*B*B.dot(E) + 1/2.0 * fPtr[RHO]*qbym2[n]/epsilon0*dt*E - 1/2.0*omega[n]*fPtr[RHO]*qbym2[n]/epsilon0*dt*B.cross(E)/2.0 )/(1.0+0.25*omega[n]*omega[n]*babs*babs);
        }
      }
      F = F + E;
      if (implicitB) {
        // this is dt*c^2*curl B
        curlBdt(0) = elfNewPtr[0] - E(0);
        curlBdt(1) = elfNewPtr[1] - E(1);
        curlBdt(2) = elfNewPtr[2] - E(2);
      }
      // if using the implicit method, add the additional curl B term to the generalised J 
      K += 1.0*curlBdt;

// fill in elements for electric field equations
      Eigen::VectorXd sol;
      if (linSolType != ANALYTIC) {
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
        // these are in the paper Smithe 2007 but are incorrect...

        /*        double coeff_e = (1.0 - 0.25*wp2 - 1.0/64.0*det2*babs*babs)/(1.0 + 0.25*wp2 + 1.0/64.0*det2*babs*babs); 
        double coeff_k = 1.0/(1.0 + 0.25*wp2 + 1.0/64.0*det2*babs*babs);
        double coeff_dot = (1.0/64.0*det2 - 1.0/16.0*gamma2)/(1.0 + 0.25*wp2 + 1.0/64.0*det2*babs*babs)/(1.0 + 0.25*wp2 + 1.0/16.0*det2*gamma2*babs*babs);
        double coeff_cross = 0.125*delta/(1.0 + 0.25*wp2 + 1.0/64.0*det2*babs*babs)/(1.0 + 0.25*wp2);     
        */
        // this is correct -- or at least I think is correct -- Jonathan Ng 2016
        double denom = ((1+theta/4)*(1+theta/4) + delta*delta*babs*babs/64.0) * (1 + theta/4 + gamma2/16*babs*babs);        
        double coeff_e = ((1+theta/4)*(1+theta/4) + gamma2*babs*babs/16*(1+theta/4));
        //        double coeff_k = 1.0/denom*0; 
        double coeff_dot = (1.0/64.0*delta*delta - 1.0/16.0*gamma2*(1+theta/4));
        double coeff_cross = 0.125*delta*(1 + gamma2*babs*babs/16 + theta/4);
        Eigen::Vector3d fk = 1.0*F + K/epsilon0; // vector holding sum of K + 2 epsilon0 E used in many calculations
        //        fk(0) = K[0] + emPtr[EX]*2*epsilon0;
        //        fk(1) = K[1] + emPtr[EY]*2*epsilon0; 
        //        fk(2) = K[2] + emPtr[EZ]*2*epsilon0; 
        emPtr[EX] = (coeff_e*fk[EX] + coeff_dot*B(0)*B.dot(fk) + coeff_cross*(B.cross(fk)(0)) )/denom; 
        emPtr[EY] = (coeff_e*fk[EY] + coeff_dot*B(1)*B.dot(fk) + coeff_cross*(B.cross(fk)(1)) )/denom; 
        emPtr[EZ] = (coeff_e*fk[EZ] + coeff_dot*B(2)*B.dot(fk) + coeff_cross*(B.cross(fk)(2)) )/denom; 

        // update the stored E field to E_n+1/2
        E(0) = (emPtr[EX] + E(0))/2.0;
        E(1) = (emPtr[EY] + E(1))/2.0;
        E(2) = (emPtr[EZ] + E(2))/2.0;

      }
      double chargeDens = 0.0;
      double keold = 0.0;
// update solution for fluids (solution is at half-time step)
      for (unsigned n=0; n<nFluids; ++n)
      {
        fluids[n]->setPtr(fPtr, idx);

// accumulate charge density from this fluid
        chargeDens += qbym[n]*fPtr[RHO];

// compute old kinetic energy before momenta are over-written
        keold = 0.5*(fPtr[RHOUX]*fPtr[RHOUX] + fPtr[RHOUY]*fPtr[RHOUY] + fPtr[RHOUZ]*fPtr[RHOUZ])/fPtr[RHO];

// momentum equation (sol has currents, so divide out charge and
// multiply by mass to give momentum density)
        if (linSolType != ANALYTIC) {
          fPtr[RHOUX] = 2*sol(fidx(n,X))/qbym[n] - fPtr[RHOUX];
          fPtr[RHOUY] = 2*sol(fidx(n,Y))/qbym[n] - fPtr[RHOUY];
          fPtr[RHOUZ] = 2*sol(fidx(n,Z))/qbym[n] - fPtr[RHOUZ];
        } else {
          // solve it analytically
          // there is another error in smithe's paper: all terms must be divided by 1/(1+.25 wc^2 dt^2), not just the first term
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
        if (hasPressure)
        {
// energy equation (there is no explicit energy source, so just
// recompute new kinetic energy to update total energy)
          fPtr[ER] = fPtr[ER] - keold
            + 0.5*(fPtr[RHOUX]*fPtr[RHOUX] + fPtr[RHOUY]*fPtr[RHOUY] + fPtr[RHOUZ]*fPtr[RHOUZ])/fPtr[RHO];
        }
      }

// update electric field
      if (linSolType != ANALYTIC){
        emPtr[EX] = 2*sol(eidx(X)) - emPtr[EX];
        emPtr[EY] = 2*sol(eidx(Y)) - emPtr[EY];
        emPtr[EZ] = 2*sol(eidx(Z)) - emPtr[EZ];
      } else {

      }
      double crhoc = chi_e*chargeDens/epsilon0;
// update electric field error potential source
      if (damp_e > 0)
        emPtr[PHIE] = (emPtr[PHIE]*damp_e-crhoc)/damp_e*std::exp(-damp_e*dt) + crhoc/damp_e;
      else
        emPtr[PHIE] += dt*crhoc;
// update magnetic field error potential source
      emPtr[PHIM] = emPtr[PHIM]*std::exp(-damp_m*dt);
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

