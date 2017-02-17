/**
 * @file	LcDGImplicitFiveMomentSrcUpdater.cpp
 *
 * @brief	Implicit updater for 5-moment source terms using a discontinuous Galerkin method
 */

// lucee includes
#include <LcDGImplicitFiveMomentSrcUpdater.h>
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

// number of components for EM (8) and fluid (5) variables
  static const unsigned numEMcomp = 8;
  static const unsigned numFldcomp = 5;

// set ids for module system
  template <> const char *DGImplicitFiveMomentSrcUpdater<1>::id = "DGImplicitFiveMomentSrc1D";
  template <> const char *DGImplicitFiveMomentSrcUpdater<2>::id = "DGImplicitFiveMomentSrc2D";
  template <> const char *DGImplicitFiveMomentSrcUpdater<3>::id = "DGImplicitFiveMomentSrc3D";

  template <unsigned NDIM>
  void
  DGImplicitFiveMomentSrcUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      basis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("DGImplicitFiveMomentSrcUpdater::readInput: Must specify configuration-space basis using 'basis'");

    nFluids = (int) tbl.getNumber("numFluids");

// charge on each fluid species
    charge = tbl.getNumVec("charge");
    if (charge.size() != nFluids)
    {
      Lucee::Except lce("DGImplicitFiveMomentSrcUpdater::readInput: 'charge' table should have exactly ");
      lce << nFluids << " entries. Provided " << charge.size() << " instead.";
      throw lce;
    }

// mass of each fluid species
    mass = tbl.getNumVec("mass");
    if (mass.size() != nFluids)
    {
      Lucee::Except lce("DGImplicitFiveMomentSrcUpdater::readInput: 'mass' table should have exactly ");
      lce << nFluids << " entries. Provided " << mass.size() << " instead.";
      throw lce;
    }

// permittivity of free space
    epsilon0 = tbl.getNumber("epsilon0");

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
  DGImplicitFiveMomentSrcUpdater<NDIM>::initialize()
  {
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  DGImplicitFiveMomentSrcUpdater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    unsigned nlocal = basis->getNumNodes();

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

    std::vector<double> zeros(6*nlocal);
    for (unsigned i=0; i<6*nlocal; ++i) zeros[i] = 0.0;
    Lucee::ConstFieldPtr<double> staticEmPtr(zeros);

    int idx[NDIM];
    Lucee::Region<NDIM, int> localRgn = emField.getRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);

      emField.setPtr(emPtr, idx);

      if (hasStatic)
        staticField->setPtr(staticEmPtr, idx);

      std::vector<double> lambda(nFluids);
      std::vector<double> omega(nFluids);

      Eigen::Vector3d B; 
      Eigen::Vector3d E;

      for (unsigned i=0; i<nlocal; ++i)
      {
        double gamma2 = 0.0;
        double delta = 0.0;
        double theta = 0.0; 

        Eigen::Vector3d F(0,0,0); 
        Eigen::Vector3d K(0,0,0);// the k vector used to update the implicit solution in Smithe(2007)

        double bx = (emPtr[i*numEMcomp+BX] + staticEmPtr[i*6+BX]);
        double by = (emPtr[i*numEMcomp+BY] + staticEmPtr[i*6+BY]);
        double bz = (emPtr[i*numEMcomp+BZ] + staticEmPtr[i*6+BZ]);

        B(0) = bx; B(1) = by; B(2) = bz;
        E(0) = emPtr[i*numEMcomp+EX]; E(1) = emPtr[i*numEMcomp+EY]; E(2) = emPtr[i*numEMcomp+EZ];
        double babs = std::sqrt(B(0)*B(0) + B(1)*B(1) + B(2)*B(2));

// fill elements corresponding to each fluid
        for (unsigned n=0; n<nFluids; ++n)
        {
          fluids[n]->setPtr(fPtr, idx);

          // use the analytic solver developed by Smithe (2007)
          lambda[n] = 1.00; //collisional parameter
          omega[n] = qbym[n]*dt; // omega/B
          double wp2 = fPtr[i*numFldcomp+RHO]*qbym2[n]/epsilon0*dt*dt;

          gamma2 += wp2*qbym2[n]*dt*dt/(1+0.25*qbym2[n]*babs*babs*dt*dt);
          delta += wp2*qbym[n]*dt/(1+0.25*qbym2[n]*babs*babs*dt*dt);
          theta += wp2/(1+qbym2[n]*dt*dt/4*babs*babs);
          Eigen::Vector3d j(fPtr[i*numFldcomp+RHOUX]*qbym[n], fPtr[i*numFldcomp+RHOUY]*qbym[n], fPtr[i*numFldcomp+RHOUZ]*qbym[n]);
          double bdotj = B.dot(j);
          Eigen::Vector3d bcrossj = B.cross(j);

          K -= dt/2.0 *(1.0 + lambda[n])*(1.0*j + 0.25*omega[n]*omega[n]*B * bdotj - 0.5*omega[n]*bcrossj )/(1.0+0.25*omega[n]*omega[n]*babs*babs);
          F -= dt/2.0*(1/4.0*omega[n]*omega[n]*fPtr[i*numFldcomp+RHO]*qbym2[n]/epsilon0*dt/2.0*B*B.dot(E) + 1/2.0 * fPtr[i*numFldcomp+RHO]*qbym2[n]/epsilon0*dt*E - 1/2.0*omega[n]*fPtr[i*numFldcomp+RHO]*qbym2[n]/epsilon0*dt*B.cross(E)/2.0 )/(1.0+0.25*omega[n]*omega[n]*babs*babs);

        }
        F = F + E;
// fill in elements for electric field equations
        Eigen::VectorXd sol;

        // this is correct -- or at least I think is correct -- Jonathan Ng 2016
        double denom = ((1+theta/4)*(1+theta/4) + delta*delta*babs*babs/64.0) * (1 + theta/4 + gamma2/16*babs*babs);        
        double coeff_e = ((1+theta/4)*(1+theta/4) + gamma2*babs*babs/16*(1+theta/4));

        double coeff_dot = (1.0/64.0*delta*delta - 1.0/16.0*gamma2*(1+theta/4));
        double coeff_cross = 0.125*delta*(1 + gamma2*babs*babs/16 + theta/4);
        Eigen::Vector3d fk = 1.0*F + K/epsilon0; // vector holding sum of K + 2 epsilon0 E used in many calculations

        emPtr[i*numEMcomp+EX] = (coeff_e*fk[EX] + coeff_dot*B(0)*B.dot(fk) + coeff_cross*(B.cross(fk)(0)) )/denom; 
        emPtr[i*numEMcomp+EY] = (coeff_e*fk[EY] + coeff_dot*B(1)*B.dot(fk) + coeff_cross*(B.cross(fk)(1)) )/denom; 
        emPtr[i*numEMcomp+EZ] = (coeff_e*fk[EZ] + coeff_dot*B(2)*B.dot(fk) + coeff_cross*(B.cross(fk)(2)) )/denom; 

        // update the stored E field to E_n+1/2
        E(0) = (emPtr[i*numEMcomp+EX] + E(0))/2.0;
        E(1) = (emPtr[i*numEMcomp+EY] + E(1))/2.0;
        E(2) = (emPtr[i*numEMcomp+EZ] + E(2))/2.0;

        double chargeDens = 0.0;
        double keold = 0.0;
// update solution for fluids (solution is at half-time step)
        for (unsigned n=0; n<nFluids; ++n)
        {
          fluids[n]->setPtr(fPtr, idx);

// accumulate charge density from this fluid
          chargeDens += qbym[n]*fPtr[i*numFldcomp+RHO];

// compute old kinetic energy before momenta are over-written
          keold = 0.5*(fPtr[i*numFldcomp+RHOUX]*fPtr[i*numFldcomp+RHOUX] + fPtr[i*numFldcomp+RHOUY]*fPtr[i*numFldcomp+RHOUY] + fPtr[i*numFldcomp+RHOUZ]*fPtr[i*numFldcomp+RHOUZ])/fPtr[i*numFldcomp+RHO];

// momentum equation (sol has currents, so divide out charge and
// multiply by mass to give momentum density)

          // solve it analytically
          // there is another error in smithe's paper: all terms must be divided by 1/(1+.25 wc^2 dt^2), not just the first term
          double j_coeff = (lambda[n] - 0.25*omega[n]*omega[n]*babs*babs)/(1.0+0.25*omega[n]*omega[n]*babs*babs);
          double f_coeff = fPtr[i*numFldcomp+RHO]*qbym2[n]/epsilon0*dt/(1.0+0.25*omega[n]*omega[n]*babs*babs);
          
          double dot_coeff = 0.25*omega[n]*omega[n]/(1.0+0.25*omega[n]*omega[n]*babs*babs);
          double cross_coeff = omega[n]/2.0/(1.0+0.25*omega[n]*omega[n]*babs*babs);
          // update E with the new solution
          Eigen::Vector3d j(fPtr[i*numFldcomp+RHOUX]*qbym[n],fPtr[i*numFldcomp+RHOUY]*qbym[n],fPtr[i*numFldcomp+RHOUZ]*qbym[n]);
          Eigen::Vector3d jf = (1.0+lambda[n])*j + E*epsilon0*dt*fPtr[i*numFldcomp+RHO]*qbym2[n]/epsilon0; // There is an error in smithe's paper here: F in equation (15) must be multiplied by a factor of wp^2 dt or the units are wrong.
          Eigen::Vector3d solj = j_coeff*j + f_coeff*E*epsilon0 + dot_coeff*B*B.dot(jf) - cross_coeff*B.cross(jf); 

          fPtr[i*numFldcomp+RHOUX] = solj(0)/qbym[n];
          fPtr[i*numFldcomp+RHOUY] = solj(1)/qbym[n];
          fPtr[i*numFldcomp+RHOUZ] = solj(2)/qbym[n];
 
          if (hasPressure)
          {
// energy equation (there is no explicit energy source, so just
// recompute new kinetic energy to update total energy)
            fPtr[i*numFldcomp+ER] = fPtr[i*numFldcomp+ER] - keold
              + 0.5*(fPtr[i*numFldcomp+RHOUX]*fPtr[i*numFldcomp+RHOUX] + fPtr[i*numFldcomp+RHOUY]*fPtr[i*numFldcomp+RHOUY] + fPtr[i*numFldcomp+RHOUZ]*fPtr[i*numFldcomp+RHOUZ])/fPtr[i*numFldcomp+RHO];
          }
        }

// update electric field

        double crhoc = chi_e*chargeDens/epsilon0;
// update electric field error potential source
        if (damp_e > 0)
          emPtr[i*numEMcomp+PHIE] = (emPtr[i*numEMcomp+PHIE]*damp_e-crhoc)/damp_e*std::exp(-damp_e*dt) + crhoc/damp_e;
        else
          emPtr[i*numEMcomp+PHIE] += dt*crhoc;
// update magnetic field error potential source
        emPtr[i*numEMcomp+PHIM] = emPtr[i*numEMcomp+PHIM]*std::exp(-damp_m*dt);
      }

    }
    
    return Lucee::UpdaterStatus();
  }
  

  template <unsigned NDIM>
  void
  DGImplicitFiveMomentSrcUpdater<NDIM>::declareTypes()
  {
    this->setLastOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class DGImplicitFiveMomentSrcUpdater<1>;
  template class DGImplicitFiveMomentSrcUpdater<2>;
  template class DGImplicitFiveMomentSrcUpdater<3>;
}
