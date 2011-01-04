/**
 * @file	LcRectCurlUpdater.cpp
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
#include <LcCurlEval.h>
#include <LcField.h>
#include <LcRectCurlUpdater.h>
#include <LcStructuredGridBase.h>

namespace Lucee
{
  template <> const char *RectCurlUpdater<1>::id = "Curl1D";
  template <> const char *RectCurlUpdater<2>::id = "Curl2D";
  template <> const char *RectCurlUpdater<3>::id = "Curl3D";

  template <unsigned NDIM>
  RectCurlUpdater<NDIM>::RectCurlUpdater()
    : UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  RectCurlUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// read multiplication factor
    alpha = tbl.getNumber("alpha");
  }

  template <unsigned NDIM>
  void
  RectCurlUpdater<NDIM>::initialize()
  {
// call base class method
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  RectCurlUpdater<NDIM>::update(double t)
  {
// get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
// get input/output arrays to compute A = B + dt*curl(V)
    const Lucee::Field<NDIM, double>& B = this->getInp<Lucee::Field<NDIM, double> >(0);
    const Lucee::Field<NDIM, double>& V = this->getInp<Lucee::Field<NDIM, double> >(1);
    Lucee::Field<NDIM, double>& A = this->getOut<Lucee::Field<NDIM, double> >(0);
    
// A <- B
    A.copy(B);

// create pointers to fields
    Lucee::ConstFieldPtr<double> Bptr = B.createConstPtr();
    Lucee::ConstFieldPtr<double> Vptr = V.createConstPtr();
    Lucee::ConstFieldPtr<double> Vptrr = V.createConstPtr();
    Lucee::FieldPtr<double> Aptr = A.createPtr();

// local region to index
    Lucee::Region<NDIM, int> localRgn = grid.getLocalBox();
// extend it to include on ghost cell on "right" of each direction
    int lg[NDIM], ug[NDIM];
    for (unsigned i=0; i<NDIM; ++i)
    {
      lg[i] = 0;
      ug[i] = 1;
    }

// time-step
    double dt = t-this->getCurrTime();
// compute factors in each direction
    double adtdx[NDIM];
    for (unsigned i=0; i<NDIM; ++i)
      adtdx[i] = alpha*dt/grid.getDx(i);

    int idx[NDIM];
// create sequencer over extended region
    Lucee::RowMajorSequencer<NDIM> seq(localRgn.extend(lg, ug));
// loop, compute curl and updating field
    while (seq.step())
    {
// get index and set pointers to appropriate places
      seq.fillWithIndex(idx);
      A.setPtr(Aptr, idx);
      B.setPtr(Bptr, idx);
      V.setPtr(Vptr, idx);
 
// compute curl in X-direction
      idx[0] += 1; // bump to one cell right
      V.setPtr(Vptrr, idx);
      idx[0] -= 1;
      Lucee::CurlEval<0>::eval(adtdx[0], &Vptr[0], &Vptrr[0], &Aptr[0]);
      if (NDIM>1)
      {
        idx[1] += 1;
        V.setPtr(Vptrr, idx);
        Lucee::CurlEval<1>::eval(adtdx[0], &Vptr[0], &Vptrr[0], &Aptr[0]);
        idx[1] -= 1;
      }
      if (NDIM>2)
      {
        idx[2] += 1;
        V.setPtr(Vptrr, idx);
        Lucee::CurlEval<2>::eval(adtdx[0], &Vptr[0], &Vptrr[0], &Aptr[0]);
        idx[2] -= 1;
      }
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  RectCurlUpdater<NDIM>::declareTypes()
  {
// two input fields
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
// one output field
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class RectCurlUpdater<1>;
  template class RectCurlUpdater<2>;
  template class RectCurlUpdater<3>;
}
