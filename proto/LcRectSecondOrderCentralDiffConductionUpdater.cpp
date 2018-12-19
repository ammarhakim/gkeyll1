/**
 * @file	LcRectSecondOrderCentralDiffConductionUpdater.cpp
 *
 * @brief	Compute 2nd order central-differences on a rectangular grid.for anisotropic heat conductivity with B field
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcRectCartGrid.h>
#include <LcRectSecondOrderCentralDiffConductionUpdater.h>
#include <LcStructGridField.h>

namespace Lucee
{
/** Class id: this is used by registration system */
  template <> const char *RectSecondOrderCentralDiffConductionUpdater<1>::id = "RectSecondOrderCentralDiffConduction1D";
  template <> const char *RectSecondOrderCentralDiffConductionUpdater<2>::id = "RectSecondOrderCentralDiffConduction2D";
  template <> const char *RectSecondOrderCentralDiffConductionUpdater<3>::id = "RectSecondOrderCentralDiffConduction3D";


  static const unsigned P11 = 0;
  static const unsigned P12 = 1;
  static const unsigned P13 = 2;
  static const unsigned P22 = 3;
  static const unsigned P23 = 4;
  static const unsigned P33 = 5;

  static const unsigned T11 = 0;
  static const unsigned T12 = 1;
  static const unsigned T13 = 2;
  static const unsigned T22 = 3;
  static const unsigned T23 = 4;
  static const unsigned T33 = 5;

  static const unsigned Q111 = 0;
  static const unsigned Q112 = 1;
  static const unsigned Q113 = 2;
  static const unsigned Q122 = 3;
  static const unsigned Q123 = 4;
  static const unsigned Q133 = 5;
  static const unsigned Q222 = 6;
  static const unsigned Q223 = 7;
  static const unsigned Q233 = 8;
  static const unsigned Q333 = 9;

  static const int BX = 3;
  static const int BY = 4;
  static const int BZ = 5;

  template <unsigned NDIM>
  RectSecondOrderCentralDiffConductionUpdater<NDIM>::RectSecondOrderCentralDiffConductionUpdater()
    : UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  RectSecondOrderCentralDiffConductionUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);
    
    alpha = tbl.getNumber("alpha"); // length scale 1/|k|
    rhoFactor = tbl.getNumber("gyroradiusFactor"); // factor to multiply 1/rho in denominator
    charge = tbl.getNumber("charge");
    mass = tbl.getNumber("mass");
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  RectSecondOrderCentralDiffConductionUpdater<NDIM>::update(double t)
  {
// get input/output fields
    const Lucee::Field<NDIM, double>& inFld = this->getInp<Lucee::Field<NDIM, double> >(0);
    const Lucee::Field<NDIM, double>& rhoFld = this->getInp<Lucee::Field<NDIM, double> >(1);
    const Lucee::Field<NDIM, double>& emField = this->getInp<Lucee::Field<NDIM, double> >(2);
    Lucee::Field<NDIM, double>& outFld = this->getOut<Lucee::Field<NDIM, double> >(0);
    Lucee::Field<NDIM, double>& qFld = this->getOut<Lucee::Field<NDIM, double> >(1);

    if (NDIM == 1)
      computeConduction1D(inFld, rhoFld, emField, outFld, qFld);
    else if (NDIM == 2)
      computeConduction2D(inFld, rhoFld, emField, outFld, qFld);
    else if (NDIM == 3)
      computeConduction3D(inFld, rhoFld, emField, outFld, qFld);

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  RectSecondOrderCentralDiffConductionUpdater<NDIM>::computeConduction1D(
                                                                       const Lucee::Field<NDIM, double>& inFld, const Lucee::Field<NDIM, double>& rhoFld, const Lucee::Field<NDIM, double>& emField, Lucee::Field<NDIM, double>& outFld, Lucee::Field<NDIM, double>& qFld)
  {
// get hold of grid
    const Lucee::RectCartGrid<1>& grid 
      = this->getGrid<Lucee::RectCartGrid<1> >();

    double dx = grid.getDx(0);
    double q[10], dT1[6], dT2[6], dT3[6], p[6], B[3], rho;
    for (unsigned i=0; i < 6; ++i){
      dT2[i] = 0; dT3[i] = 0;
    }
    
    // first, calculate q
    for (int i=inFld.getLower(0)-1; i<inFld.getUpper(0)+1; ++i) {
      // we need a single ghost cell for q
      B[0] = emField(i,BX);
      B[1] = emField(i,BY);
      B[2] = emField(i,BZ);
      rho = rhoFld(i,0);
      
      for (unsigned k=0; k<inFld.getNumComponents(); ++k){
        p[k] = inFld(i,k);
        dT1[k] = 0.5*(inFld(i+1,k)/rhoFld(i,0) - inFld(i-1,k)/rhoFld(i-1,0))/dx; 
      }
      heatFlux(p, dT1, dT2, dT3, B, rho, q);
      for (unsigned m = 0; m < qFld.getNumComponents(); ++m) {
        qFld(i,m) = q[m];
      }
    }

    // Now calculate Divergence of Q 
    // d_l q_ij l = d_1 q_ij1 + d_2 q_ij 2 + d_3 q_ij3
    for (int i=inFld.getLower(0); i<inFld.getUpper(0); ++i) {
      outFld(i,P11) = 0.5*(qFld(i+1,Q111) - qFld(i-1,Q111))/dx;
      outFld(i,P12) = 0.5*(qFld(i+1,Q112) - qFld(i-1,Q112))/dx;
      outFld(i,P13) = 0.5*(qFld(i+1,Q113) - qFld(i-1,Q113))/dx;
      outFld(i,P22) = 0.5*(qFld(i+1,Q122) - qFld(i-1,Q122))/dx;
      outFld(i,P23) = 0.5*(qFld(i+1,Q123) - qFld(i-1,Q123))/dx;
      outFld(i,P33) = 0.5*(qFld(i+1,Q133) - qFld(i-1,Q133))/dx;
    }
  }

  template <unsigned NDIM>
  void
  RectSecondOrderCentralDiffConductionUpdater<NDIM>::computeConduction2D(
                                                                       const Lucee::Field<NDIM, double>& inFld, const Lucee::Field<NDIM, double>& rhoFld, const Lucee::Field<NDIM, double>& emField, Lucee::Field<NDIM, double>& outFld, Lucee::Field<NDIM, double>& qFld)
  {
// get hold of grid
    const Lucee::RectCartGrid<2>& grid 
      = this->getGrid<Lucee::RectCartGrid<2> >();

    double dx = grid.getDx(0);
    double dy = grid.getDx(1);

    double q[10], dT1[6], dT2[6], dT3[6], p[6], B[3], rho;
    for (unsigned i=0; i < 6; ++i){
      dT3[i] = 0;
    }

    // first, calculate q
    for (int i=inFld.getLower(0)-1; i<inFld.getUpper(0)+1; ++i) {
      for (int j=inFld.getLower(1)-1; j<inFld.getUpper(1)+1; ++j) {
      // we need a single ghost cell for q
        B[0] = emField(i,j,BX);
        B[1] = emField(i,j,BY);
        B[2] = emField(i,j,BZ);
        rho = rhoFld(i,j,0);
        
        for (unsigned k=0; k<inFld.getNumComponents(); ++k){
          p[k] = inFld(i,j,k);
          dT1[k] = 0.5*(inFld(i+1,j,k)/rhoFld(i+1,j,0) - inFld(i-1,j,k)/rhoFld(i-1,j,0))/dx;
          dT2[k] = 0.5*(inFld(i,j+1,k)/rhoFld(i,j+1,0) - inFld(i,j-1,k)/rhoFld(i,j-1,0))/dy;
          /*          dT1[k] = (inFld(i,j,k)/rhoFld(i,j,0) - inFld(i-1,j,k)/rhoFld(i-1,j,0))/dx;
                      dT2[k] = (inFld(i,j,k)/rhoFld(i,j,0) - inFld(i,j-1,k)/rhoFld(i,j-1,0))/dy;*/
        }
        heatFlux(p, dT1, dT2, dT3, B, rho, q);
        for (unsigned m = 0; m < qFld.getNumComponents(); ++m) {
          qFld(i,j,m) = q[m];
        }
      }
    }

    for (int i=inFld.getLower(0); i<inFld.getUpper(0); ++i) {
      for (int j=inFld.getLower(1); j<inFld.getUpper(1); ++j) {
        outFld(i,j,P11) = 0.5*(qFld(i+1,j,Q111) - qFld(i-1,j,Q111))/dx + 0.5*(qFld(i,j+1,Q112) - qFld(i,j-1,Q112))/dy;
        outFld(i,j,P12) = 0.5*(qFld(i+1,j,Q112) - qFld(i-1,j,Q112))/dx + 0.5*(qFld(i,j+1,Q122) - qFld(i,j-1,Q122))/dy;
        outFld(i,j,P13) = 0.5*(qFld(i+1,j,Q113) - qFld(i-1,j,Q113))/dx + 0.5*(qFld(i,j+1,Q123) - qFld(i,j-1,Q123))/dy;
        outFld(i,j,P22) = 0.5*(qFld(i+1,j,Q122) - qFld(i-1,j,Q122))/dx + 0.5*(qFld(i,j+1,Q222) - qFld(i,j-1,Q222))/dy;
        outFld(i,j,P23) = 0.5*(qFld(i+1,j,Q123) - qFld(i-1,j,Q123))/dx + 0.5*(qFld(i,j+1,Q223) - qFld(i,j-1,Q223))/dy;
        outFld(i,j,P33) = 0.5*(qFld(i+1,j,Q133) - qFld(i-1,j,Q133))/dx + 0.5*(qFld(i,j+1,Q233) - qFld(i,j-1,Q233))/dy;

      }
    }

  }

  template <unsigned NDIM>
  void
  RectSecondOrderCentralDiffConductionUpdater<NDIM>::computeConduction3D(
                                                                       const Lucee::Field<NDIM, double>& inFld, const Lucee::Field<NDIM, double>& rhoFld, const Lucee::Field<NDIM, double>& emField, Lucee::Field<NDIM, double>& outFld, Lucee::Field<NDIM, double>& qFld)
  {
// get hold of grid
    const Lucee::RectCartGrid<3>& grid 
      = this->getGrid<Lucee::RectCartGrid<3> >();

    double dx = grid.getDx(0);
    double dy = grid.getDx(1);
    double dz = grid.getDx(2);


    double q[10], dT1[6], dT2[6], dT3[6], p[6], B[3], rho;

    // first, calculate q
    for (int i=inFld.getLower(0)-1; i<inFld.getUpper(0)+1; ++i) {
      for (int j=inFld.getLower(1)-1; j<inFld.getUpper(1)+1; ++j) {
        for (int l=inFld.getLower(2)-1; l<inFld.getUpper(2)+1; ++l) {
          B[0] = emField(i,j,l,BX);
          B[1] = emField(i,j,l,BY);
          B[2] = emField(i,j,l,BZ);
          rho = rhoFld(i,j,l,0);
      // we need a single ghost cell for q
          for (unsigned k=0; k<inFld.getNumComponents(); ++k){
            //            std::cout<<"hello world " << i << " " << l <<" \n";
            p[k] = inFld(i,j,l,k);
            dT1[k] = 0.5*(inFld(i+1,j,l,k)/rhoFld(i+1,j,l,0) - inFld(i-1,j,l,k)/rhoFld(i-1,j,l,0))/dx;
            dT2[k] = 0.5*(inFld(i,j+1,l,k)/rhoFld(i,j+1,l,0) - inFld(i,j-1,l,k)/rhoFld(i,j-1,l,0))/dy;
            dT3[k] = 0.5*(inFld(i,j,l+1,k)/rhoFld(i,j,l+1,0) - inFld(i,j,l-1,k)/rhoFld(i,j,l-1,0))/dz;
          }
          heatFlux(p, dT1, dT2, dT3, B, rho, q);
          for (unsigned m = 0; m < qFld.getNumComponents(); ++m) {
            qFld(i,j,l,m) = q[m];
          }
        }
      }
    }

    for (int i=inFld.getLower(0)-1; i<inFld.getUpper(0)+1; ++i) {
      for (int j=inFld.getLower(1)-1; j<inFld.getUpper(1)+1; ++j) {
        for (int l=inFld.getLower(2)-1; l<inFld.getUpper(2)+1; ++l) {
          outFld(i,j,l,P11) = 0.5*((qFld(i+1,j,l,Q111) - qFld(i-1,j,l,Q111))/dx + (qFld(i,j+1,l,Q112) - qFld(i,j-1,l,Q112))/dy + (qFld(i,j,l+1,Q113) - qFld(i,j,l-1,Q113))/dz);
          outFld(i,j,l,P12) = 0.5*((qFld(i+1,j,l,Q112) - qFld(i-1,j,l,Q112))/dx + (qFld(i,j+1,l,Q122) - qFld(i,j-1,l,Q122))/dy + (qFld(i,j,l+1,Q123) - qFld(i,j,l-1,Q123))/dz);
          outFld(i,j,l,P13) = 0.5*((qFld(i+1,j,l,Q113) - qFld(i-1,j,l,Q113))/dx + (qFld(i,j+1,l,Q123) - qFld(i,j-1,l,Q123))/dy + (qFld(i,j,l+1,Q133) - qFld(i,j,l-1,Q133))/dz);
          outFld(i,j,l,P22) = 0.5*((qFld(i+1,j,l,Q122) - qFld(i-1,j,l,Q122))/dx + (qFld(i,j+1,l,Q222) - qFld(i,j-1,l,Q222))/dy + (qFld(i,j,l+1,Q223) - qFld(i,j,l-1,Q223))/dz);
          outFld(i,j,l,P23) = 0.5*((qFld(i+1,j,l,Q123) - qFld(i-1,j,l,Q123))/dx + (qFld(i,j+1,l,Q223) - qFld(i,j-1,l,Q223))/dy + (qFld(i,j,l+1,Q233) - qFld(i,j,l-1,Q233))/dz);
          outFld(i,j,l,P33) = 0.5*((qFld(i+1,j,l,Q133) - qFld(i-1,j,l,Q133))/dx + (qFld(i,j+1,l,Q233) - qFld(i,j-1,l,Q233))/dy + (qFld(i,j,l+1,Q333) - qFld(i,j,l-1,Q333))/dz);
        }
      }
    }


  }



  template <unsigned NDIM>
  void
  RectSecondOrderCentralDiffConductionUpdater<NDIM>::declareTypes()
  {
    this->setLastInpVarType(typeid(Lucee::Field<NDIM, double>));
    //    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->setLastOutVarType(typeid(Lucee::Field<NDIM, double>));
  }


  template <unsigned NDIM>
  void 
  RectSecondOrderCentralDiffConductionUpdater<NDIM>::heatFlux(const double* p, const double* dT1, const double* dT2, const double* dT3, const double* B, const double rho, double* q) { 
    // this is the form of the collisional heat flux in the regularised moment method
    // we are using q = qpar + qperp*factor
    // q_isotropic = d_i T jk * rho * chi, symmetrised
    // q_parallel = d_parallel_i T_jk
    // q_perp = q_isotropic - q_parallel
    // ideally this will go to the isotropic limit as B -> 0

    // The hammett perkins heat flux is n vt sqrt(8/pi) partial_[k T_ij] where v_t is sqrt(P/rho)

    double Babs = std::sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]) + 1e-8;
    double b1 = B[0]/Babs; double b2 = B[1]/Babs; double b3 = B[2]/Babs;
    double vt = std::sqrt(2*(p[P11]+p[P22]+p[P33])/rho);
    // use the inverse larmor radius to avoid divide by zero
    double invLarmorRadius = std::abs(charge)*Babs/mass/vt;
    
    double chipar = alpha*vt;
    // scale chipar by 1/(1+1/kperp^2rho^2)
    double chiperp = chipar/(1+invLarmorRadius*invLarmorRadius*alpha*alpha*rhoFactor); 
    // instead of alpha, one would want to use the magnetic field length sclae. 

    std::vector<double> qpar(10);
    std::vector<double> qiso(10);
    // parallel heat flux, automatically generated
    qpar[Q111] = rho*(b1*(b1*dT1[T11] + b2*dT2[T11] + b3*dT3[T11]) + b1*(b1*dT1[T11] + b2*dT2[T11] + b3*dT3[T11]) + b1*(b1*dT1[T11] + b2*dT2[T11] + b3*dT3[T11]));
    qpar[Q112] = rho*(b1*(b1*dT1[T12] + b2*dT2[T12] + b3*dT3[T12]) + b1*(b1*dT1[T12] + b2*dT2[T12] + b3*dT3[T12]) + b2*(b1*dT1[T11] + b2*dT2[T11] + b3*dT3[T11]));
    qpar[Q113] = rho*(b1*(b1*dT1[T13] + b2*dT2[T13] + b3*dT3[T13]) + b1*(b1*dT1[T13] + b2*dT2[T13] + b3*dT3[T13]) + b3*(b1*dT1[T11] + b2*dT2[T11] + b3*dT3[T11]));
    qpar[Q122] = rho*(b1*(b1*dT1[T22] + b2*dT2[T22] + b3*dT3[T22]) + b2*(b1*dT1[T12] + b2*dT2[T12] + b3*dT3[T12]) + b2*(b1*dT1[T12] + b2*dT2[T12] + b3*dT3[T12]));
    qpar[Q123] = rho*(b1*(b1*dT1[T23] + b2*dT2[T23] + b3*dT3[T23]) + b2*(b1*dT1[T13] + b2*dT2[T13] + b3*dT3[T13]) + b3*(b1*dT1[T12] + b2*dT2[T12] + b3*dT3[T12]));
    qpar[Q133] = rho*(b1*(b1*dT1[T33] + b2*dT2[T33] + b3*dT3[T33]) + b3*(b1*dT1[T13] + b2*dT2[T13] + b3*dT3[T13]) + b3*(b1*dT1[T13] + b2*dT2[T13] + b3*dT3[T13]));
    qpar[Q222] = rho*(b2*(b1*dT1[T22] + b2*dT2[T22] + b3*dT3[T22]) + b2*(b1*dT1[T22] + b2*dT2[T22] + b3*dT3[T22]) + b2*(b1*dT1[T22] + b2*dT2[T22] + b3*dT3[T22]));
    qpar[Q223] = rho*(b2*(b1*dT1[T23] + b2*dT2[T23] + b3*dT3[T23]) + b2*(b1*dT1[T23] + b2*dT2[T23] + b3*dT3[T23]) + b3*(b1*dT1[T22] + b2*dT2[T22] + b3*dT3[T22]));
    qpar[Q233] = rho*(b2*(b1*dT1[T33] + b2*dT2[T33] + b3*dT3[T33]) + b3*(b1*dT1[T23] + b2*dT2[T23] + b3*dT3[T23]) + b3*(b1*dT1[T23] + b2*dT2[T23] + b3*dT3[T23]));
    qpar[Q333] = rho*(b3*(b1*dT1[T33] + b2*dT2[T33] + b3*dT3[T33]) + b3*(b1*dT1[T33] + b2*dT2[T33] + b3*dT3[T33]) + b3*(b1*dT1[T33] + b2*dT2[T33] + b3*dT3[T33]));
    // isotropic heat flux, automatically generated
    qiso[Q111] = rho*(dT1[T11] + dT1[T11] + dT1[T11]);
    qiso[Q112] = rho*(dT1[T12] + dT1[T12] + dT2[T11]);
    qiso[Q113] = rho*(dT1[T13] + dT1[T13] + dT3[T11]);
    qiso[Q122] = rho*(dT1[T22] + dT2[T12] + dT2[T12]);
    qiso[Q123] = rho*(dT1[T23] + dT2[T13] + dT3[T12]);
    qiso[Q133] = rho*(dT1[T33] + dT3[T13] + dT3[T13]);
    qiso[Q222] = rho*(dT2[T22] + dT2[T22] + dT2[T22]);
    qiso[Q223] = rho*(dT2[T23] + dT2[T23] + dT3[T22]);
    qiso[Q233] = rho*(dT2[T33] + dT3[T23] + dT3[T23]);
    qiso[Q333] = rho*(dT3[T33] + dT3[T33] + dT3[T33]);
    // perp heat flux, automatically generated. 
    q[Q111] = (chipar*qpar[Q111] + chiperp*(qiso[Q111] - qpar[Q111]))/3.0;
    q[Q112] = (chipar*qpar[Q112] + chiperp*(qiso[Q112] - qpar[Q112]))/3.0;
    q[Q113] = (chipar*qpar[Q113] + chiperp*(qiso[Q113] - qpar[Q113]))/3.0;
    q[Q122] = (chipar*qpar[Q122] + chiperp*(qiso[Q122] - qpar[Q122]))/3.0;
    q[Q123] = (chipar*qpar[Q123] + chiperp*(qiso[Q123] - qpar[Q123]))/3.0;
    q[Q133] = (chipar*qpar[Q133] + chiperp*(qiso[Q133] - qpar[Q133]))/3.0;
    q[Q222] = (chipar*qpar[Q222] + chiperp*(qiso[Q222] - qpar[Q222]))/3.0;
    q[Q223] = (chipar*qpar[Q223] + chiperp*(qiso[Q223] - qpar[Q223]))/3.0;
    q[Q233] = (chipar*qpar[Q233] + chiperp*(qiso[Q233] - qpar[Q233]))/3.0;
    q[Q333] = (chipar*qpar[Q333] + chiperp*(qiso[Q333] - qpar[Q333]))/3.0;
  }


// instantiations
  template class RectSecondOrderCentralDiffConductionUpdater<1>;
  template class RectSecondOrderCentralDiffConductionUpdater<2>;
  template class RectSecondOrderCentralDiffConductionUpdater<3>;
}
