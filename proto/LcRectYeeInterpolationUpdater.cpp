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
#include <LcRectCartGrid.h>
#include <LcRectYeeInterpolationUpdater.h>
#include <LcStructGridField.h>

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
      // It is generally necessary to do this as the last FV cell will 
      // use an edge/face outside the regular domain
      std::vector<double> gup = tbl.getNumVec("ghostUpdate");
      if (gup.size() != 2)
        throw Lucee::Except(
          "RectYeeInterpolationUpdater:readInput: The 'ghostUpdate' table should have exactly two numbers");
      ghostUpdates[0] = (unsigned) gup[0];
      ghostUpdates[1] = (unsigned) gup[1];
    }
    else
    {
// by default do not update any ghost region
      ghostUpdates[0] = ghostUpdates[1] = 0; 
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
    const Lucee::RectCartGrid<1>& grid 
      = this->getGrid<Lucee::RectCartGrid<1> >();

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
          cdFld(i,0) = 0.5*(inFld(i,0) + inFld(i-1,0));
          cdFld(i,1) = inFld(i,1);
          cdFld(i,2) = inFld(i,2);
          cdFld(i,3) = inFld(i,3);
          cdFld(i,4) = 0.5*( inFld(i,4) + inFld(i-1,4) );
          cdFld(i,5) = 0.5*( inFld(i,5) + inFld(i-1,5) );
        }
      }
    } else {
      // Regular Yee cell with the Electric fields on the edges and the 
      // magnetic field on the cell centres.
      if (fwd) {// interpolate to cell centre
        for (int i=inFld.getLower(0)+ghostUpdates[0]; i<inFld.getUpper(0)+ghostUpdates[1]; ++i) {
          cdFld(i,3) = 0.5*(inFld(i,3) + inFld(i+1,3));
          cdFld(i,4) = inFld(i,4);
          cdFld(i,5) = inFld(i,5);
          cdFld(i,0) = inFld(i,0);
          cdFld(i,1) = 0.5*( inFld(i,1) + inFld(i+1,1) );
          cdFld(i,2) = 0.5*( inFld(i,2) + inFld(i+1,2) );
        }
      } else { // interpolate back to edge/face
        for (int i=inFld.getLower(0)+ghostUpdates[0]; i<inFld.getUpper(0)+ghostUpdates[1]; ++i){
          cdFld(i,3) = 0.5*(inFld(i,3) + inFld(i-1,3));
          cdFld(i,4) = inFld(i,4);
          cdFld(i,5) = inFld(i,5);
          cdFld(i,0) = inFld(i,0);
          cdFld(i,1) = 0.5*( inFld(i,1) + inFld(i-1,1) );
          cdFld(i,2) = 0.5*( inFld(i,2) + inFld(i-1,2) );
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
    const Lucee::RectCartGrid<2>& grid 
      = this->getGrid<Lucee::RectCartGrid<2> >();

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
            
            cdFld(i,j,0) = 0.5*(inFld(i,j,0) + inFld(i-1,j,0));
            cdFld(i,j,1) = 0.5*(inFld(i,j,1) + inFld(i,j-1,1));
            cdFld(i,j,2) = inFld(i,j,2);
            cdFld(i,j,3) = 0.5*( inFld(i,j,3) + inFld(i,j-1,3) );
            cdFld(i,j,4) = 0.5*( inFld(i,j,4) + inFld(i-1,j,4) );
            cdFld(i,j,5) = 0.25*( inFld(i,j,5) + inFld(i-1,j,5) + inFld(i,j-1,5) + inFld(i-1,j-1,5) );
          }
      }
    } else { // Regular Yee cell
      if (fwd) {// interpolate to cell centre
        for (int i=inFld.getLower(0)+ghostUpdates[0]; i<inFld.getUpper(0)+ghostUpdates[1]; ++i)
          for (int j=inFld.getLower(1)+ghostUpdates[0]; j<inFld.getUpper(1)+ghostUpdates[1]; ++j){
            // in the dual yee cell, electric fields are at faces and 
            // magnetic fields are at edges
            cdFld(i,j,3) = 0.5*(inFld(i,j,3) + inFld(i+1,j,3));
            cdFld(i,j,4) = 0.5*(inFld(i,j,4) + inFld(i,j+1,4));
            cdFld(i,j,5) = inFld(i,j,5);
            cdFld(i,j,0) = 0.5*( inFld(i,j,0) + inFld(i,j+1,0) );
            cdFld(i,j,1) = 0.5*( inFld(i,j,1) + inFld(i+1,j,1) );
            cdFld(i,j,2) = 0.25*( inFld(i,j,2) + inFld(i+1,j,2) + inFld(i,j+1,2) + inFld(i+1,j+1,2) );
            
          }
      } else { // interpolate back to edge/face
        for (int i=inFld.getLower(0)+ghostUpdates[0]; i<inFld.getUpper(0)+ghostUpdates[1]; ++i)
          for (int j=inFld.getLower(1)+ghostUpdates[0]; j<inFld.getUpper(1)+ghostUpdates[1]; ++j){
            
            cdFld(i,j,3) = 0.5*(inFld(i,j,3) + inFld(i-1,j,3));
            cdFld(i,j,4) = 0.5*(inFld(i,j,4) + inFld(i,j-1,4));
            cdFld(i,j,5) = inFld(i,j,5);
            cdFld(i,j,0) = 0.5*( inFld(i,j,0) + inFld(i,j-1,0) );
            cdFld(i,j,1) = 0.5*( inFld(i,j,1) + inFld(i-1,j,1) );
            cdFld(i,j,2) = 0.25*( inFld(i,j,2) + inFld(i-1,j,2) + inFld(i,j-1,2) + inFld(i-1,j-1,2) );
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
    const Lucee::RectCartGrid<3>& grid 
      = this->getGrid<Lucee::RectCartGrid<3> >();

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
              
              cdFld(i,j,l,0) = 0.5*(inFld(i,j,l,0) + inFld(i-1,j,l,0));
              cdFld(i,j,l,1) = 0.5*(inFld(i,j,l,1) + inFld(i,j-1,l,1));
              cdFld(i,j,l,2) = 0.5*(inFld(i,j,l,2) + inFld(i,j,l-1,2));
              cdFld(i,j,l,3) = 0.25*( inFld(i,j,l,3) + inFld(i,j-1,l,3) + inFld(i,j,l-1,3) + inFld(i,j-1,l-1,3) );
              cdFld(i,j,l,4) = 0.25*( inFld(i,j,l,4) + inFld(i-1,j,l,4) + inFld(i,j,l-1,4) + inFld(i-1,j,l-1,4) );
              cdFld(i,j,l,5) = 0.25*( inFld(i,j,l,5) + inFld(i-1,j,l,5) + inFld(i,j-1,l,5) + inFld(i-1,j-1,l,5) );
            }
      }
    } else { // Regular  Yee cell
      if (fwd) {// interpolate to cell centre
        for (int i=inFld.getLower(0)+ghostUpdates[0]; i<inFld.getUpper(0)+ghostUpdates[1]; ++i)
          for (int j=inFld.getLower(1)+ghostUpdates[0]; j<inFld.getUpper(1)+ghostUpdates[1]; ++j)
            for (int l=inFld.getLower(2)+ghostUpdates[0]; l<inFld.getUpper(2)+ghostUpdates[1]; ++l){
          
            // in the dual yee cell, electric fields are at faces and 
            // magnetic fields are at edges
              cdFld(i,j,l,3) = 0.5*(inFld(i,j,l,3) + inFld(i+1,j,l,3));
              cdFld(i,j,l,4) = 0.5*(inFld(i,j,l,4) + inFld(i,j+1,l,4));
              cdFld(i,j,l,5) = 0.5*(inFld(i,j,l,5) + inFld(i,j,l+1,5));
              cdFld(i,j,l,0) = 0.25*( inFld(i,j,l,0) + inFld(i,j+1,l,0) + inFld(i,j,l+1,0) + inFld(i,j+1,l+1,0) );
              cdFld(i,j,l,1) = 0.25*( inFld(i,j,l,1) + inFld(i+1,j,l,1) + inFld(i,j,l+1,1) + inFld(i+1,j,l+1,1) );
              cdFld(i,j,l,2) = 0.25*( inFld(i,j,l,2) + inFld(i+1,j,l,2) + inFld(i,j+1,l,2) + inFld(i+1,j+1,l,2) );
              
            }
      } else { // interpolate back to edge/face
        for (int i=inFld.getLower(0)+ghostUpdates[0]; i<inFld.getUpper(0)+ghostUpdates[1]; ++i)
          for (int j=inFld.getLower(1)+ghostUpdates[0]; j<inFld.getUpper(1)+ghostUpdates[1]; ++j)
            for (int l=inFld.getLower(2)+ghostUpdates[0]; l<inFld.getUpper(2)+ghostUpdates[1]; ++l){
              
              cdFld(i,j,l,3) = 0.5*(inFld(i,j,l,3) + inFld(i-1,j,l,3));
              cdFld(i,j,l,4) = 0.5*(inFld(i,j,l,4) + inFld(i,j-1,l,4));
              cdFld(i,j,l,5) = 0.5*(inFld(i,j,l,5) + inFld(i,j,l-1,5));
              cdFld(i,j,l,0) = 0.25*( inFld(i,j,l,0) + inFld(i,j-1,l,0) + inFld(i,j,l-1,0) + inFld(i,j-1,l-1,0) );
              cdFld(i,j,l,1) = 0.25*( inFld(i,j,l,1) + inFld(i-1,j,l,1) + inFld(i,j,l-1,1) + inFld(i-1,j,l-1,1) );
              cdFld(i,j,l,2) = 0.25*( inFld(i,j,l,2) + inFld(i-1,j,l,2) + inFld(i,j-1,l,2) + inFld(i-1,j-1,l,2) );
            }
      }
    }
  }

// instantiations
  template class RectYeeInterpolationUpdater<1>;
  template class RectYeeInterpolationUpdater<2>;
  template class RectYeeInterpolationUpdater<3>;
}
