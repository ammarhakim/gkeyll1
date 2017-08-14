/**
 * @file	LcNonUniEdgeFaceCurlUpdater.cpp
 *
 * @brief	Compute curl on rectangular grids.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcNonUniEdgeFaceCurlUpdater.h>
#include <LcStructuredGridBase.h>

namespace Lucee
{
  template <> const char *NonUniEdgeFaceCurlUpdater<1>::id = "NonUniEdgeFaceCurl1D";
  template <> const char *NonUniEdgeFaceCurlUpdater<2>::id = "NonUniEdgeFaceCurl2D";
  template <> const char *NonUniEdgeFaceCurlUpdater<3>::id = "NonUniEdgeFaceCurl3D";

  template <unsigned NDIM>
  NonUniEdgeFaceCurlUpdater<NDIM>::NonUniEdgeFaceCurlUpdater()
    : UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  NonUniEdgeFaceCurlUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);
// read stuff needed to run updater
    alpha = tbl.getNumber("alpha");
    speed = tbl.getNumber("speed");
    speed = speed > 0 ? speed : -speed;
    cfl = tbl.getNumber("cfl");
    if (tbl.hasNumVec("ghostUpdate"))
    {
      std::vector<double> gup = tbl.getNumVec("ghostUpdate");
      if (gup.size() != 2)
        throw Lucee::Except(
          "NonUniEdgeFaceCurlUpdater:readInput: The 'ghostUpdate' table should have exactly two numbers");
      ghostUpdates[0] = (unsigned) gup[0];
      ghostUpdates[1] = (unsigned) gup[1];
    }
    else
    {
// by default do not update any ghost region
      ghostUpdates[0] = ghostUpdates[1] = 0; 
    }

// get hold of grid
// for the non uniform grid we want to find the min dx in each direction
// this is so that we can set the cfl condition correctly 
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    int idx[NDIM];
    // array to hold the left and right vertices
    double vl[NDIM]; double vr[NDIM]; 
    for (unsigned dim = 0; dim < NDIM; ++dim){
      // initialise the index to a valid number
      idx[dim] = localRgn.getLower(dim);
    }
    for (unsigned dim = 0; dim < NDIM; ++dim) {
      dxMin[dim] = -1.0;
      int sliceLower = localRgn.getLower(dim);
      int sliceUpper = localRgn.getUpper(dim);
      for (unsigned j = sliceLower; j < sliceUpper; ++j){
        idx[dim] = j;
        grid.setIndex(idx);
        // get right vertex
        grid.getVertex(vr);
        idx[dim] = j-1;
        // get left vertex
        // there is no need to reset the index after each direction
        // because this is a stretched cartesian grid
        grid.setIndex(idx);
        grid.getVertex(vl);
        double dx = std::abs(vr[dim]-vl[dim]);
        if (dxMin[dim] > dx || dxMin[dim] < 0) {
          dxMin[dim] = dx;
        }
      }
    }


  }

  template <unsigned NDIM>
  void
  NonUniEdgeFaceCurlUpdater<NDIM>::initialize()
  {
// call base class method
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  NonUniEdgeFaceCurlUpdater<NDIM>::update(double t)
  {
// get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

// time-step
    double dt = t-this->getCurrTime();

    double cfla = 0.0;
// first ensure that CFL condition is met
    for (unsigned d=0; d<NDIM; ++d)
      cfla = std::max(cfla, speed*std::abs(dt)/dxMin[d]);

    if (cfla > 1.01*cfl) 
    { // time-step too large (the 1.01 is used to avoid thrashing around dt)
      double newDt = dt*cfl/cfla;
      return UpdaterStatus(false, newDt);
    }

// arrays to store the vertex coordinates
    double vl[NDIM]; double vr[NDIM];

// get input/output arrays to compute A = B + dt*alpha*curl(V). A and
// B must be edge-centered and V must be face centered.
    const Lucee::Field<NDIM, double>& B = this->getInp<Lucee::Field<NDIM, double> >(0);
    const Lucee::Field<NDIM, double>& V = this->getInp<Lucee::Field<NDIM, double> >(1);
    Lucee::Field<NDIM, double>& A = this->getOut<Lucee::Field<NDIM, double> >(0);

// create pointers to fields
    Lucee::ConstFieldPtr<double> Vptr = V.createConstPtr();
    Lucee::ConstFieldPtr<double> Vptrl = V.createConstPtr();
    Lucee::FieldPtr<double> Aptr = A.createPtr();
    
// local region to index
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    Lucee::Region<NDIM, int> extRgn = A.getExtRegion();

// A <- B
    A.copy(B);
// loop, updating slices in each direction
    for (unsigned dir=0; dir<NDIM; ++dir)
    {

// create sequencer to loop over *each* 1D slice in 'dir' direction
      Lucee::RowMajorSequencer<NDIM> seq(extRgn.deflate(dir));
//      Lucee::RowMajorSequencer<NDIM> seqExt(extRgn.deflate(dir));
// lower and upper bounds of 1D slice
      int sliceLower = localRgn.getLower(dir) + ghostUpdates[0];
      int sliceUpper = localRgn.getUpper(dir) + ghostUpdates[1];

      double adtdx, dxLocal; 
// loop over each 1D slice
      while (seq.step())
      {
        int idx[NDIM], idxl[NDIM];
        seq.fillWithIndex(idx);
        seq.fillWithIndex(idxl);
// (if-test is outside loop over slice to amortize its cost)
        if (dir == 0)
        { // X-direction update
          for (int i=sliceLower; i<sliceUpper; ++i)
          {
            idx[dir] = i; // current cell
            idxl[dir] = i-1; // left cell
// set pointers to proper location          
            A.setPtr(Aptr, idx);
            V.setPtr(Vptr, idx);
            V.setPtr(Vptrl, idxl);

// in this updater we calculate dx by using the centroids of the left
// and current cell
            // get right vertex
            grid.setIndex(idx);
            grid.getCentroid(vr);
            // get left vertex
            grid.setIndex(idxl);
            grid.getCentroid(vl);
            dxLocal = vr[dir] - vl[dir];
            adtdx = alpha*dt/dxLocal;

            Aptr[1] += -adtdx*(Vptr[2]-Vptrl[2]);
            Aptr[2] += adtdx*(Vptr[1]-Vptrl[1]);
          }
        }
        else if (dir == 1)
        { // Y-direction update
          for (int i=sliceLower; i<sliceUpper; ++i)
          {
            idx[dir] = i; // current cell
            idxl[dir] = i-1; // left cell
// set pointers to proper location          
            A.setPtr(Aptr, idx);
            V.setPtr(Vptr, idx);
            V.setPtr(Vptrl, idxl);

// in this updater we calculate dx by using the centroids of the left
// and current cell
            // get right vertex
            grid.setIndex(idx);
            grid.getCentroid(vr);
            // get left vertex
            grid.setIndex(idxl);
            grid.getCentroid(vl);
            dxLocal = vr[dir] - vl[dir];
            adtdx = alpha*dt/dxLocal;

            Aptr[0] += adtdx*(Vptr[2]-Vptrl[2]);
            Aptr[2] += -adtdx*(Vptr[0]-Vptrl[0]);
          }
        }
        else
        { // Z-direction update
          for (int i=sliceLower; i<sliceUpper; ++i)
          {
            idx[dir] = i; // current cell
            idxl[dir] = i-1; // left cell
// set pointers to proper location          
            A.setPtr(Aptr, idx);
            V.setPtr(Vptr, idx);
            V.setPtr(Vptrl, idxl);

// in this updater we calculate dx by using the centroids of the left
// and current cell
            // get right vertex
            grid.setIndex(idx);
            grid.getCentroid(vr);
            // get left vertex
            grid.setIndex(idxl);
            grid.getCentroid(vl);
            dxLocal = vr[dir] - vl[dir];
            adtdx = alpha*dt/dxLocal;

            Aptr[0] += -adtdx*(Vptr[1]-Vptrl[1]);
            Aptr[1] += adtdx*(Vptr[0]-Vptrl[0]);
          }
        }
      }
    }

    return Lucee::UpdaterStatus(true, dt*cfl/cfla);
  }

  template <unsigned NDIM>
  void
  NonUniEdgeFaceCurlUpdater<NDIM>::declareTypes()
  {
// two input fields
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
// one output field
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class NonUniEdgeFaceCurlUpdater<1>;
  template class NonUniEdgeFaceCurlUpdater<2>;
  template class NonUniEdgeFaceCurlUpdater<3>;
}
