/**
 * @file	LcRectYeeInterpolationUpdater.cpp
 *
 * @brief	Interpolate fields from regular FV cells to Yee cells and vice versa
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcRectYeeInterpolationUpdater.h>
#include <LcStructGridField.h>
#include <LcStructuredGridBase.h>

namespace Lucee
{
/** Class id: this is used by registration system */
  template <> const char *RectYeeInterpolationUpdater<1>::id = "RectYeeInterpolation1D";
  template <> const char *RectYeeInterpolationUpdater<2>::id = "RectYeeInterpolation2D";
  template <> const char *RectYeeInterpolationUpdater<3>::id = "RectYeeInterpolation3D";

  template <unsigned NDIM>
  RectYeeInterpolationUpdater<NDIM>::RectYeeInterpolationUpdater()
    : UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  RectYeeInterpolationUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);
    if (tbl.hasBool("toCentre")){
      fwd = tbl.getBool("toCentre");
    } else {
        Lucee::Except lce("RectYeeInterpolationUpdater::readInput: 'toCentre' must be true or false");
        throw lce;
    }
    if (tbl.hasBool("isDual")){
      dual = tbl.getBool("isDual");
    } else {
        Lucee::Except lce("RectYeeInterpolationUpdater::readInput: 'isDual  ' must be true or false");
        throw lce;

    }

    if (tbl.hasNumVec("ghostUpdate"))
    {
      // This is sometimes important when the final edge cell must be updated
      // and is not set by an explicit boundary condition. 
      std::vector<double> gup = tbl.getNumVec("ghostUpdate");
      if (gup.size() != 2)
        throw Lucee::Except(
          "RectYeeInterpolationUpdater:readInput: The 'ghostUpdate' table should have exactly two numbers");
      ghostUpdates[0] = gup[0];
      ghostUpdates[1] = gup[1];
    }
    else
    {
// by default do not update any ghost region
      ghostUpdates[0] = ghostUpdates[1] = 0; 
    }

    // 
    // When using perfectly conducting boundaries, the interpolation of E can 
    // cause a force imbalance at the boundary. This leads to an inward
    // propagating grid-scale oscillation. 
    // This effect is also present in the Dual cell but doesn't appear
    // to be as obvious (likely because the reflected quantities do not
    // live on the boundary) 
    // This fix is only implemented for the regular cell and is off by default
    // Set bcL/bcU [i] to zero if not intending to apply to said boundary. 
   

    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    Lucee::Region<NDIM, int> globalRgn = grid.getGlobalRegion();

    for (unsigned i = 0; i < NDIM; ++i) {
      // search for boundaries. 
      if (localRgn.getLower(i) == globalRgn.getLower(i)) {
        globalLower[i] = globalRgn.getLower(i);
        // this is a lower boundary
      }
      if (localRgn.getUpper(i) == globalRgn.getUpper(i)) {
        globalUpper[i] = globalRgn.getUpper(i); 
        // this is an upper boundary 
      }
    }

    if (tbl.hasNumVec("copyLowerSkinCell"))
    {
      std::vector<double> lowerSkin = tbl.getNumVec("copyLowerSkinCell");
      if (lowerSkin.size() != NDIM)
        throw Lucee::Except(
          "RectYeeInterpolationUpdater:readInput: The 'lowerSkin' table should have exactly NDIM numbers");
      for (int i = 0; i < NDIM; ++i)
        bcL[i] = lowerSkin[i];
      
    } else {
      for (int i = 0; i < NDIM; ++i) {
        // by default do not copy the last cell
        bcL[i] = 0;
      }
    }
    if (tbl.hasNumVec("copyUpperSkinCell"))
    {
      std::vector<double> upperSkin = tbl.getNumVec("copyUpperSkinCell");
      if (upperSkin.size() != NDIM)
        throw Lucee::Except(
          "RectYeeInterpolationUpdater:readInput: The 'lowerSkin' table should have exactly NDIM numbers");
      for (int i = 0; i < NDIM; ++i)
        bcU[i] = upperSkin[i];
    } else {
      for (int i = 0; i < NDIM; ++i) { 
        // by default do not copy the last cell
        bcU[i] = 0;
      }
    }
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  RectYeeInterpolationUpdater<NDIM>::update(double t)
  {
// get input/output fields
    const Lucee::Field<NDIM, double>& inFld = this->getInp<Lucee::Field<NDIM, double> >(0);
    Lucee::Field<NDIM, double>& cdFld = this->getOut<Lucee::Field<NDIM, double> >(0);

    if (NDIM == 1)
      computeCentralDifference1D(inFld, cdFld);
    else if (NDIM == 2)
      computeCentralDifference2D(inFld, cdFld);
    else if (NDIM == 3)
      computeCentralDifference3D(inFld, cdFld);

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  RectYeeInterpolationUpdater<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void
  RectYeeInterpolationUpdater<NDIM>::computeCentralDifference1D(
    const Lucee::Field<NDIM, double>& inFld, Lucee::Field<NDIM, double>& cdFld)
  {
// get hold of grid
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();

    // when going to cell centre there is no need to find dx
    // only when going back to the edges do we need the distances
    // for the bilinear interpolation

    // we are working in at most a 2d plane. 
    // The general formula in 1-d is 
    // f(x) = (x2-x)/dx f(x1) + (x-x1)/dx f(x2)
    int idx[NDIM];
    double xc[1], xl[1], xr[1]; 
    
    if (dual) {
      if (fwd) {// interpolate to cell centre
        for (int i=inFld.getLower(0)+ghostUpdates[0]; i<inFld.getUpper(0)+ghostUpdates[1]; ++i) {
          // in the dual yee cell, electric fields are at faces and 
          // magnetic fields are at edges
          cdFld(i,0) = 0.5*(inFld(i,0) + inFld(i+1,0));
          cdFld(i,1) = inFld(i,1);
          cdFld(i,2) = inFld(i,2);
          cdFld(i,3) = inFld(i,3);
          cdFld(i,4) = 0.5*( inFld(i,4) + inFld(i+1,4) );
          cdFld(i,5) = 0.5*( inFld(i,5) + inFld(i+1,5) );
        }
      } else { // interpolate back to edge/face
        for (int i=inFld.getLower(0)+ghostUpdates[0]; i<inFld.getUpper(0)+ghostUpdates[1]; ++i){
          idx[0] = i;
          grid.setIndex(idx);
          grid.getCentroid(xr);
          grid.getVertex(xc);

          idx[0] = i-1;
          grid.setIndex(idx);
          grid.getCentroid(xl);
          //assert( xc[0] > xl[0] && xc[0]  <xr[0]);
          double dx = xr[0] - xl[0];
          double dxl = (xc[0] - xl[0])/dx;
          double dxr = (xr[0] - xc[0])/dx;
          cdFld(i,0) = (inFld(i,0)*dxl + inFld(i-1,0)*dxr);
          cdFld(i,1) = inFld(i,1);
          cdFld(i,2) = inFld(i,2);
          cdFld(i,3) = inFld(i,3);
          cdFld(i,4) = ( inFld(i,4)*dxl + inFld(i-1,4)*dxr );
          cdFld(i,5) = ( inFld(i,5)*dxl + inFld(i-1,5)*dxr );
        }
      }
    } else {
      // Regular Yee cell with the Electric fields on the edges and the 
      // magnetic field on the cell centres.
      if (fwd) {// interpolate to cell centre
        for (int i=inFld.getLower(0)+ghostUpdates[0]; i<inFld.getUpper(0)+ghostUpdates[1]; ++i) {

          int xLower = (i == globalLower[0])*bcL[0];
          int xUpper = (i == globalUpper[0]-1)*bcU[0];
          
          int x0 = (1-xLower)*(1+xUpper); int x1 = (1+xLower)*(1-xUpper);
            
          cdFld(i,3) = 0.5*(inFld(i,3)*x0 + inFld(i+1,3)*x1);
          cdFld(i,4) = inFld(i,4);
          cdFld(i,5) = inFld(i,5);
          cdFld(i,0) = inFld(i,0);
          cdFld(i,1) = 0.5*( inFld(i,1)*x0 + inFld(i+1,1)*x1 );
          cdFld(i,2) = 0.5*( inFld(i,2)*x0 + inFld(i+1,2)*x1 );
        }
      } else { // interpolate back to edge/face
        for (int i=inFld.getLower(0)+ghostUpdates[0]; i<inFld.getUpper(0)+ghostUpdates[1]; ++i){

          idx[0] = i;
          grid.setIndex(idx);
          grid.getCentroid(xr);
          grid.getVertex(xc);

          idx[0] = i-1;
          grid.setIndex(idx);
          grid.getCentroid(xl);
          //assert( xc[0] > xl[0] && xc[0]  <xr[0]);
          double dx = xr[0] - xl[0];
          double dxl = (xc[0] - xl[0])/dx;
          double dxr = (xr[0] - xc[0])/dx;

          cdFld(i,3) = (inFld(i,3)*dxl + inFld(i-1,3)*dxr);
          cdFld(i,4) = inFld(i,4);
          cdFld(i,5) = inFld(i,5);
          cdFld(i,0) = inFld(i,0);
          cdFld(i,1) = ( inFld(i,1)*dxl + inFld(i-1,1)*dxr );
          cdFld(i,2) = ( inFld(i,2)*dxl + inFld(i-1,2)*dxr );
        }
      }

    }
  }



  template <unsigned NDIM>
  void
  RectYeeInterpolationUpdater<NDIM>::computeCentralDifference2D(
    const Lucee::Field<NDIM, double>& inFld, Lucee::Field<NDIM, double>& cdFld)
  {
// get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();

    // we are working in at most a 2d plane. 
    // The general formula in 1-d is 
    // f(x) = (x2-x)/dx f(x1) + (x-x1)/dx f(x2)
    // The general formula in 2d is 
    // f(x,y) = (f(x1,y1)*dxr*dyr + f(x2,y1)*dxl*dyr + f(x1,y2)*dxr*dyl + f(x2,y2)*dxl*dyl)/dx/dy
    int idx[NDIM];
    double xc[2], xl[2], xr[2]; 

    if (dual) {
      if (fwd) {// interpolate to cell centre
        for (int i=inFld.getLower(0)+ghostUpdates[0]; i<inFld.getUpper(0)+ghostUpdates[1]; ++i)
          for (int j=inFld.getLower(1)+ghostUpdates[0]; j<inFld.getUpper(1)+ghostUpdates[1]; ++j){
            // in the dual yee cell, electric fields are at faces and 
            // magnetic fields are at edges
            cdFld(i,j,0) = 0.5*(inFld(i,j,0) + inFld(i+1,j,0));
            cdFld(i,j,1) = 0.5*(inFld(i,j,1) + inFld(i,j+1,1));
            cdFld(i,j,2) = inFld(i,j,2);
            cdFld(i,j,3) = 0.5*( inFld(i,j,3) + inFld(i,j+1,3) );
            cdFld(i,j,4) = 0.5*( inFld(i,j,4) + inFld(i+1,j,4) );
            cdFld(i,j,5) = 0.25*( inFld(i,j,5) + inFld(i+1,j,5) + inFld(i,j+1,5) + inFld(i+1,j+1,5) );
            
          }
      } else { // interpolate back to edge/face
        for (int i=inFld.getLower(0)+ghostUpdates[0]; i<inFld.getUpper(0)+ghostUpdates[1]; ++i)
          for (int j=inFld.getLower(1)+ghostUpdates[0]; j<inFld.getUpper(1)+ghostUpdates[1]; ++j){
            // in 2-d we need to calculate dx, dy, dxl, dyl, dxr, dyr
            idx[0] = i; idx[1] = j;
            grid.setIndex(idx);
            grid.getCentroid(xr);
            grid.getVertex(xc);

            idx[0] = i-1; idx[1] = j-1;
            grid.setIndex(idx);
            grid.getCentroid(xl);
            //assert( xc[0] > xl[0] && xc[0]  <xr[0]);
            double dx = xr[0] - xl[0];
            double dxl = (xc[0] - xl[0])/dx;
            double dxr = (xr[0] - xc[0])/dx;
            double dy = xr[1] - xl[1];
            double dyl = (xc[1] - xl[1])/dy;
            double dyr = (xr[1] - xc[1])/dy;

            cdFld(i,j,0) = (inFld(i,j,0)*dxl + inFld(i-1,j,0)*dxr);
            cdFld(i,j,1) = (inFld(i,j,1)*dyl + inFld(i,j-1,1)*dyr);
            cdFld(i,j,2) = inFld(i,j,2);
            cdFld(i,j,3) = ( inFld(i,j,3)*dyl + inFld(i,j-1,3)*dyr );
            cdFld(i,j,4) = ( inFld(i,j,4)*dxl + inFld(i-1,j,4)*dxr );
            cdFld(i,j,5) = ( inFld(i,j,5)*dxl*dyl + inFld(i-1,j,5)*dxr*dyl + inFld(i,j-1,5)*dxl*dyr + inFld(i-1,j-1,5)*dxr*dyr );
          }
      }
    } else { // Regular Yee cell
      // Regular Yee cell with the electric fields on the edges and the 
      // magnetic field on the cell centres.

      if (fwd) {// interpolate to cell centre
        for (int i=inFld.getLower(0)+ghostUpdates[0]; i<inFld.getUpper(0)+ghostUpdates[1]; ++i)
          for (int j=inFld.getLower(1)+ghostUpdates[0]; j<inFld.getUpper(1)+ghostUpdates[1]; ++j){

            
            int xLower = (i == globalLower[0])*bcL[0];
            int xUpper = (i == globalUpper[0]-1)*bcU[0];
            int yLower = (j == globalLower[1])*bcL[1];
            int yUpper = (j == globalUpper[1]-1)*bcU[1]; 
            
            int x0 = (1-xLower)*(1+xUpper); int x1 = (1+xLower)*(1-xUpper);
            int y0 = (1-yLower)*(1+yUpper); int y1 = (1+yLower)*(1-yUpper);

            cdFld(i,j,3) = 0.5*(inFld(i,j,3)*x0 + inFld(i+1,j,3)*x1);
            cdFld(i,j,4) = 0.5*(inFld(i,j,4)*y0 + inFld(i,j+1,4)*y1);
            cdFld(i,j,5) = inFld(i,j,5);
            cdFld(i,j,0) = 0.5*( inFld(i,j,0)*y0 + inFld(i,j+1,0)*y1 );
            cdFld(i,j,1) = 0.5*( inFld(i,j,1)*x0 + inFld(i+1,j,1)*x1 );
            cdFld(i,j,2) = 0.25*( inFld(i,j,2)*x0*y0 + inFld(i+1,j,2)*x1*y0 + inFld(i,j+1,2)*x0*y1 + inFld(i+1,j+1,2)*x1*y1 );
            
          }
      } else { // interpolate back to edge/face
        for (int i=inFld.getLower(0)+ghostUpdates[0]; i<inFld.getUpper(0)+ghostUpdates[1]; ++i)
          for (int j=inFld.getLower(1)+ghostUpdates[0]; j<inFld.getUpper(1)+ghostUpdates[1]; ++j){
            // in 2-d we need to calculate dx, dy, dxl, dyl, dxr, dyr
            idx[0] = i; idx[1] = j;
            grid.setIndex(idx);
            grid.getCentroid(xr);
            grid.getVertex(xc);

            idx[0] = i-1; idx[1] = j-1;
            grid.setIndex(idx);
            grid.getCentroid(xl);
            //assert( xc[0] > xl[0] && xc[0]  <xr[0]);
            double dx = xr[0] - xl[0];
            double dxl = (xc[0] - xl[0])/dx;
            double dxr = (xr[0] - xc[0])/dx;
            double dy = xr[1] - xl[1];
            double dyl = (xc[1] - xl[1])/dy;
            double dyr = (xr[1] - xc[1])/dy;         

            cdFld(i,j,3) = (inFld(i,j,3)*dxl + inFld(i-1,j,3)*dxr);
            cdFld(i,j,4) = (inFld(i,j,4)*dyl + inFld(i,j-1,4)*dyr);
            cdFld(i,j,5) = inFld(i,j,5);
            cdFld(i,j,0) = ( inFld(i,j,0)*dyl + inFld(i,j-1,0)*dyr );
            cdFld(i,j,1) = ( inFld(i,j,1)*dxl + inFld(i-1,j,1)*dxr );
            cdFld(i,j,2) = ( inFld(i,j,2)*dxl*dyl + inFld(i-1,j,2)*dxr*dyl + inFld(i,j-1,2)*dxl*dyr + inFld(i-1,j-1,2)*dxr*dyr );
          }
      }
    }
  }

  template <unsigned NDIM>
  void
  RectYeeInterpolationUpdater<NDIM>::computeCentralDifference3D(
    const Lucee::Field<NDIM, double>& inFld, Lucee::Field<NDIM, double>& cdFld)
  {
// get hold of grid
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    int idx[NDIM];
    double xc[3], xl[3], xr[3]; 

    if (dual) {
      if (fwd) {// interpolate to cell centre
        for (int i=inFld.getLower(0)+ghostUpdates[0]; i<inFld.getUpper(0)+ghostUpdates[1]; ++i)
          for (int j=inFld.getLower(1)+ghostUpdates[0]; j<inFld.getUpper(1)+ghostUpdates[1]; ++j)
            for (int l=inFld.getLower(2)+ghostUpdates[0]; l<inFld.getUpper(2)+ghostUpdates[1]; ++l){
          
            // in the dual yee cell, electric fields are at faces and 
            // magnetic fields are at edges
              cdFld(i,j,l,0) = 0.5*(inFld(i,j,l,0) + inFld(i+1,j,l,0));
              cdFld(i,j,l,1) = 0.5*(inFld(i,j,l,1) + inFld(i,j+1,l,1));
              cdFld(i,j,l,2) = 0.5*(inFld(i,j,l,2) + inFld(i,j,l+1,2));
              cdFld(i,j,l,3) = 0.25*( inFld(i,j,l,3) + inFld(i,j+1,l,3) + inFld(i,j,l+1,3) + inFld(i,j+1,l+1,3) );
              cdFld(i,j,l,4) = 0.25*( inFld(i,j,l,4) + inFld(i+1,j,l,4) + inFld(i,j,l+1,4) + inFld(i+1,j,l+1,4) );
              cdFld(i,j,l,5) = 0.25*( inFld(i,j,l,5) + inFld(i+1,j,l,5) + inFld(i,j+1,l,5) + inFld(i+1,j+1,l,5) );
              
            }
      } else { // interpolate back to edge/face
        for (int i=inFld.getLower(0)+ghostUpdates[0]; i<inFld.getUpper(0)+ghostUpdates[1]; ++i)
          for (int j=inFld.getLower(1)+ghostUpdates[0]; j<inFld.getUpper(1)+ghostUpdates[1]; ++j)
            for (int l=inFld.getLower(2)+ghostUpdates[0]; l<inFld.getUpper(2)+ghostUpdates[1]; ++l){

              idx[0] = i; idx[1] = j; idx[2] = l;
              grid.setIndex(idx);
              grid.getCentroid(xr);
              grid.getVertex(xc);

              idx[0] = i-1; idx[1] = j-1; idx[2] = l-1;
              grid.setIndex(idx);
              grid.getCentroid(xl);
              //assert( xc[0] > xl[0] && xc[0]  <xr[0]);
              double dx = xr[0] - xl[0];
              double dxl = (xc[0] - xl[0])/dx;
              double dxr = (xr[0] - xc[0])/dx;
              double dy = xr[1] - xl[1];
              double dyl = (xc[1] - xl[1])/dy;
              double dyr = (xr[1] - xc[1])/dy;  
              double dz = xr[2] - xl[2];
              double dzl = (xc[2] - xl[2])/dz;
              double dzr = (xr[2] - xc[2])/dz;

              cdFld(i,j,l,0) = (inFld(i,j,l,0)*dxl + inFld(i-1,j,l,0)*dxr);
              cdFld(i,j,l,1) = (inFld(i,j,l,1)*dyl + inFld(i,j-1,l,1)*dyr);
              cdFld(i,j,l,2) = (inFld(i,j,l,2)*dzl + inFld(i,j,l-1,2)*dzr);
              cdFld(i,j,l,3) = ( inFld(i,j,l,3)*dzl*dyl + inFld(i,j-1,l,3)*dzl*dyr + inFld(i,j,l-1,3)*dzr*dyl + inFld(i,j-1,l-1,3)*dzr*dyr );
              cdFld(i,j,l,4) = ( inFld(i,j,l,4)*dxl*dzl + inFld(i-1,j,l,4)*dxr*dzl + inFld(i,j,l-1,4)*dxl*dzr + inFld(i-1,j,l-1,4)*dxr*dzr );
              cdFld(i,j,l,5) = ( inFld(i,j,l,5)*dxl*dyl + inFld(i-1,j,l,5)*dxr*dyl + inFld(i,j-1,l,5)*dxl*dyr + inFld(i-1,j-1,l,5)*dxr*dyr );
            }
      }
    } else { // Regular  Yee cell
      // Regular Yee cell with the electric fields on the edges and the 
      // magnetic field on the cell centres.


      if (fwd) {// interpolate to cell centre
        for (int i=inFld.getLower(0)+ghostUpdates[0]; i<inFld.getUpper(0)+ghostUpdates[1]; ++i)
          for (int j=inFld.getLower(1)+ghostUpdates[0]; j<inFld.getUpper(1)+ghostUpdates[1]; ++j)
            for (int l=inFld.getLower(2)+ghostUpdates[0]; l<inFld.getUpper(2)+ghostUpdates[1]; ++l){
              // This logic applies the interpolation to only the first/last cells.
              int xLower = (i == globalLower[0])*bcL[0];
              int xUpper = (i == globalUpper[0]-1)*bcU[0];
              int yLower = (j == globalLower[1])*bcL[1];
              int yUpper = (j == globalUpper[1]-1)*bcU[1]; 
              int zLower = (l == globalLower[2])*bcL[2];
              int zUpper = (l == globalUpper[2]-1)*bcU[2];

              int x0 = (1-xLower)*(1+xUpper); int x1 = (1+xLower)*(1-xUpper);
              int y0 = (1-yLower)*(1+yUpper); int y1 = (1+yLower)*(1-yUpper);
              int z0 = (1-zLower)*(1+zUpper); int z1 = (1+zLower)*(1-zUpper);  

              cdFld(i,j,l,3) = 0.5*(inFld(i,j,l,3)*x0 + inFld(i+1,j,l,3)*x1);
              cdFld(i,j,l,4) = 0.5*(inFld(i,j,l,4)*y0 + inFld(i,j+1,l,4)*y1);
              cdFld(i,j,l,5) = 0.5*(inFld(i,j,l,5)*z0 + inFld(i,j,l+1,5)*z1);
              cdFld(i,j,l,0) = 0.25*( inFld(i,j,l,0)*y0*z0 + inFld(i,j+1,l,0)*y1*z0 + inFld(i,j,l+1,0)*y0*z1 + inFld(i,j+1,l+1,0)*y1*z1 );
              cdFld(i,j,l,1) = 0.25*( inFld(i,j,l,1)*x0*z0 + inFld(i+1,j,l,1)*x1*z0 + inFld(i,j,l+1,1)*x0*z1 + inFld(i+1,j,l+1,1)*x1*z1 );
              cdFld(i,j,l,2) = 0.25*( inFld(i,j,l,2)*x0*y0 + inFld(i+1,j,l,2)*x1*y0 + inFld(i,j+1,l,2)*x0*y1 + inFld(i+1,j+1,l,2)*x1*y1 );
              
            }
      } else { // interpolate back to edge/face
        for (int i=inFld.getLower(0)+ghostUpdates[0]; i<inFld.getUpper(0)+ghostUpdates[1]; ++i)
          for (int j=inFld.getLower(1)+ghostUpdates[0]; j<inFld.getUpper(1)+ghostUpdates[1]; ++j)
            for (int l=inFld.getLower(2)+ghostUpdates[0]; l<inFld.getUpper(2)+ghostUpdates[1]; ++l){
              
              idx[0] = i; idx[1] = j; idx[2] = l;
              grid.setIndex(idx);
              grid.getCentroid(xr);
              grid.getVertex(xc);

              idx[0] = i-1; idx[1] = j-1; idx[2] = l-1;
              grid.setIndex(idx);
              grid.getCentroid(xl);
              //assert( xc[0] > xl[0] && xc[0]  <xr[0]);

              double dx = xr[0] - xl[0];
              double dxl = (xc[0] - xl[0])/dx;
              double dxr = (xr[0] - xc[0])/dx;
              double dy = xr[1] - xl[1];
              double dyl = (xc[1] - xl[1])/dy;
              double dyr = (xr[1] - xc[1])/dy;  
              double dz = xr[2] - xl[2];
              double dzl = (xc[2] - xl[2])/dz;
              double dzr = (xr[2] - xc[2])/dz;

              cdFld(i,j,l,3) = (inFld(i,j,l,3)*dxl + inFld(i-1,j,l,3)*dxr);
              cdFld(i,j,l,4) = (inFld(i,j,l,4)*dyl + inFld(i,j-1,l,4)*dyr);
              cdFld(i,j,l,5) = (inFld(i,j,l,5)*dzl + inFld(i,j,l-1,5)*dzr);
              cdFld(i,j,l,0) = ( inFld(i,j,l,0)*dyl*dzl + inFld(i,j-1,l,0)*dyr*dzl + inFld(i,j,l-1,0)*dzr*dyl + inFld(i,j-1,l-1,0)*dyr*dzr );
              cdFld(i,j,l,1) = ( inFld(i,j,l,1)*dxl*dzl + inFld(i-1,j,l,1)*dxr*dzl + inFld(i,j,l-1,1)*dxl*dzr + inFld(i-1,j,l-1,1)*dxr*dzr );
              cdFld(i,j,l,2) = ( inFld(i,j,l,2)*dxl*dyl + inFld(i-1,j,l,2)*dxr*dyl + inFld(i,j-1,l,2)*dxl*dyr + inFld(i-1,j-1,l,2)*dxr*dyr );
            }
      }
    }

  }

// instantiations
  template class RectYeeInterpolationUpdater<1>;
  template class RectYeeInterpolationUpdater<2>;
  template class RectYeeInterpolationUpdater<3>;
}
