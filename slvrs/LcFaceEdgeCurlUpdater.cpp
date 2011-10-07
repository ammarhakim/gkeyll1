/**
 * @file	LcFaceEdgeCurlUpdater.cpp
 *
 * @brief	Compute curl on rectangular grids.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcFaceEdgeCurlUpdater.h>
#include <LcStructuredGridBase.h>

namespace Lucee
{
  template <> const char *FaceEdgeCurlUpdater<1>::id = "FaceEdgeCurl1D";
  template <> const char *FaceEdgeCurlUpdater<2>::id = "FaceEdgeCurl2D";
  template <> const char *FaceEdgeCurlUpdater<3>::id = "FaceEdgeCurl3D";

  template <unsigned NDIM>
  FaceEdgeCurlUpdater<NDIM>::FaceEdgeCurlUpdater()
    : UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  FaceEdgeCurlUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
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
          "EdgeFaceCurlUpdater:readInput: The 'ghostUpdate' table should have exactly two numbers");
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
  void
  FaceEdgeCurlUpdater<NDIM>::initialize()
  {
// call base class method
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  FaceEdgeCurlUpdater<NDIM>::update(double t)
  {
// get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

// time-step
    double dt = t-this->getCurrTime();

    double cfla = 0.0;
// first ensure that CFL condition is met
    for (unsigned d=0; d<NDIM; ++d)
      cfla = std::max(cfla, speed*dt/grid.getDx(d));

    if (cfla > 1.01*cfl) 
    { // time-step too large (the 1.01 is used to avoid thrashing around dt)
      double newDt = dt*cfl/cfla;
      return UpdaterStatus(false, newDt);
    }

// get input/output arrays to compute A = B + dt*alpha*curl(V). A and
// B must be face-centered and V must be edge centered.
    const Lucee::Field<NDIM, double>& B = this->getInp<Lucee::Field<NDIM, double> >(0);
    const Lucee::Field<NDIM, double>& V = this->getInp<Lucee::Field<NDIM, double> >(1);
    Lucee::Field<NDIM, double>& A = this->getOut<Lucee::Field<NDIM, double> >(0);

// create pointers to fields
    Lucee::ConstFieldPtr<double> Vptr = V.createConstPtr();
    Lucee::ConstFieldPtr<double> Vptrr = V.createConstPtr();
    Lucee::FieldPtr<double> Aptr = A.createPtr();
    
// local region to index
    Lucee::Region<NDIM, int> localRgn = grid.getLocalBox();

// A <- B
    A.copy(B);
// loop, updating slices in each direction
    for (unsigned dir=0; dir<NDIM; ++dir)
    {
      double adtdx = alpha*dt/grid.getDx(dir);
// create sequencer to loop over *each* 1D slice in 'dir' direction
      Lucee::RowMajorSequencer<NDIM> seq(localRgn.deflate(dir));

// lower and upper bounds of 1D slice
      int sliceLower = localRgn.getLower(dir) + ghostUpdates[0];
      int sliceUpper = localRgn.getUpper(dir) + ghostUpdates[1];

// loop over each 1D slice
      while (seq.step())
      {
        int idx[NDIM], idxr[NDIM];
        seq.fillWithIndex(idx);
        seq.fillWithIndex(idxr);
// (if-test is outside loop over slice to amortize its cost)
        if (dir == 0)
        { // X-direction update
          for (int i=sliceLower; i<sliceUpper; ++i)
          {
            idx[dir] = i; // current cell
            idxr[dir] = i+1; // right cell
// set pointers to proper location          
            A.setPtr(Aptr, idx);
            V.setPtr(Vptr, idx);
            V.setPtr(Vptrr, idxr);
            Aptr[1] += -adtdx*(Vptrr[2]-Vptr[2]);
            Aptr[2] += adtdx*(Vptrr[1]-Vptr[1]);
          }
        }
        else if (dir == 1)
        { // Y-direction update
          for (int i=sliceLower; i<sliceUpper; ++i)
          {
            idx[dir] = i; // current cell
            idxr[dir] = i+1; // right cell
// set pointers to proper location          
            A.setPtr(Aptr, idx);
            V.setPtr(Vptr, idx);
            V.setPtr(Vptrr, idxr);
            Aptr[0] += adtdx*(Vptrr[2]-Vptr[2]);
            Aptr[2] += -adtdx*(Vptrr[0]-Vptr[0]);
          }
        }
        else
        { // Z-direction update
          for (int i=sliceLower; i<sliceUpper; ++i)
          {
            idx[dir] = i; // current cell
            idxr[dir] = i+1; // right cell
// set pointers to proper location          
            A.setPtr(Aptr, idx);
            V.setPtr(Vptr, idx);
            V.setPtr(Vptrr, idxr);
            Aptr[0] += -adtdx*(Vptrr[1]-Vptr[1]);
            Aptr[1] += adtdx*(Vptrr[0]-Vptr[0]);
          }
        }
      }
    }

    return Lucee::UpdaterStatus(true, dt*cfl/cfla);
  }

  template <unsigned NDIM>
  void
  FaceEdgeCurlUpdater<NDIM>::declareTypes()
  {
// two input fields
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
// one output field
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class FaceEdgeCurlUpdater<1>;
  template class FaceEdgeCurlUpdater<2>;
  template class FaceEdgeCurlUpdater<3>;
}
