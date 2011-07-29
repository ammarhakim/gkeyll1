/**
 * @file   LcGridOdePointIntegrator.cpp
 *
 * @brief   ODE integrator over complete grid.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcConstFieldPtr.h>
#include <LcFieldPtr.h>
#include <LcGridOdePointIntegrator.h>

namespace Lucee
{
  template <unsigned NDIM>
  GridOdePointIntegrator<NDIM>::GridOdePointIntegrator(const Lucee::StructuredGridBase<NDIM>& grid)
    : GridOdeIntegrator<NDIM>(grid, 0)
  {
  }

  template <unsigned NDIM>
  GridOdePointIntegrator<NDIM>::~GridOdePointIntegrator()
  {
    rhs.clear();
  }

  template <unsigned NDIM>
  void
  GridOdePointIntegrator<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method first
    Lucee::GridOdeIntegrator<NDIM>::readInput(tbl);

// get list of  terms to use on RHS
    if ( !tbl.hasTable("terms") )
      Lucee::Except lce("GridOdePointIntegrator::readInput: Must specify terms to to use!");

    Lucee::LuaTable trmTbl = tbl.getTable("terms");
    rhs = trmTbl.template getAllObjects<Lucee::PointSourceIfc>();
  }

  template <unsigned NDIM>
  void
  GridOdePointIntegrator<NDIM>::integrate(double t0, double t1, Lucee::Field<NDIM, double>& sol)
  {
    ts.resize(sol.getNumComponents());
// time-step to use
    double dt = t1-t0;
// update using RK4 scheme
    rk4(dt, sol);
  }

  template <unsigned NDIM>
  void
  GridOdePointIntegrator<NDIM>::rk4(double dt, Lucee::Field<NDIM, double>& sol)
  {
// number of components to update
    unsigned n = sol.getNumComponents();

    std::vector<double> ql(n), src(n), srct(n), srcm(n);
// get reference to grid
    const Lucee::StructuredGridBase<NDIM>& grid = this->getGrid();

// create pointers to fields
    Lucee::FieldPtr<double> solPtr = sol.createPtr();
    Lucee::ConstFieldPtr<double> inpPtr = sol.createConstPtr();

    double hh = dt/2.0;
    double h6 = dt/6.0;

    int idx[NDIM];
    double xc[3];
    Lucee::RowMajorSequencer<NDIM> seq(sol.getRegion());
// loop over field, updating solution in each cell
    while (seq.step())
    {
// get index and get centroid coordinate
      seq.fillWithIndex(idx);
      grid.getCentriod(xc);
// set pointers
      sol.setPtr(inpPtr, idx);

// RK stage 1
      calcSource(xc, &inpPtr[0], src);
      for (unsigned i=0; i<n; ++i)
        ql[i] = inpPtr[i] + hh*src[i];

// RK stage 2
      calcSource(xc, &ql[0], srct);
      for (unsigned i=0; i<n; ++i)
        ql[i] = inpPtr[i] + hh*srct[i];

// RK stage 3
      calcSource(xc, &ql[0], srcm);
      for (unsigned i=0; i<n; ++i)
      {
        ql[i] = inpPtr[i] + dt*srct[i];
        srcm[i] = srct[i] + srcm[i];
      }

// RK stage 4
      calcSource(xc, &ql[0], srct);
// perform final update
      sol.setPtr(solPtr, idx);
      for (unsigned i=0; i<n; ++i)
        solPtr[i] += h6*(src[i] + srct[i] + 2*srcm[i]);
    }
  }

  template <unsigned NDIM>
  void
  GridOdePointIntegrator<NDIM>::calcSource(const double xc[3], const double *inp, std::vector<double>& src)
  {
    unsigned n = src.size();
// zap sources first
    for (unsigned k=0; k<n; ++k)
      src[k] = 0.0;
// compute each source, accumulating it
    for (unsigned i=0; i<rhs.size(); ++i)
    {
      for (unsigned k=0; k<n; ++k) 
        ts[k] = 0.0;
      rhs[i]->calcSource(xc, inp, &ts[0]);
      for (unsigned k=0; k<n; ++k)
        src[k] += ts[k];
    }
  }

// instantiations
  template class GridOdePointIntegrator<1>;
  template class GridOdePointIntegrator<2>;
  template class GridOdePointIntegrator<3>;
}
