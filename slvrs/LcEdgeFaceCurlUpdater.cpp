/**
 * @file	LcEdgeFaceCurlUpdater.cpp
 *
 * @brief	Compute curl on rectangular grids.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcEdgeFaceCurlUpdater.h>
#include <LcStructuredGridBase.h>

namespace Lucee
{
  template <> const char *EdgeFaceCurlUpdater<1>::id = "EdgeFaceCurl1D";
  template <> const char *EdgeFaceCurlUpdater<2>::id = "EdgeFaceCurl2D";
  template <> const char *EdgeFaceCurlUpdater<3>::id = "EdgeFaceCurl3D";

  template <unsigned NDIM>
  EdgeFaceCurlUpdater<NDIM>::EdgeFaceCurlUpdater()
    : UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  EdgeFaceCurlUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);
// read multiplication factor
    alpha = tbl.getNumber("alpha");
  }

  template <unsigned NDIM>
  void
  EdgeFaceCurlUpdater<NDIM>::initialize()
  {
// call base class method
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  EdgeFaceCurlUpdater<NDIM>::update(double t)
  {
// get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
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
    Lucee::Region<NDIM, int> localRgn = grid.getLocalBox();
// time-step
    double dt = t-this->getCurrTime();

// A <- B
    A.copy(B);
// loop, updating slices in each direction
    for (unsigned dir=0; dir<NDIM; ++dir)
    {
      double adtdx = alpha*dt/grid.getDx(dir);
// create sequencer to loop over *each* 1D slice in 'dir' direction
      Lucee::RowMajorSequencer<NDIM> seq(localRgn.deflate(dir));

// lower and upper bounds of 1D slice
      int sliceLower = localRgn.getLower(dir);
      int sliceUpper = localRgn.getUpper(dir)+1;

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
            Aptr[0] += -adtdx*(Vptr[1]-Vptrl[1]);
            Aptr[1] += adtdx*(Vptr[0]-Vptrl[0]);
          }
        }
      }
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  EdgeFaceCurlUpdater<NDIM>::declareTypes()
  {
// two input fields
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
// one output field
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class EdgeFaceCurlUpdater<1>;
  template class EdgeFaceCurlUpdater<2>;
  template class EdgeFaceCurlUpdater<3>;
}
