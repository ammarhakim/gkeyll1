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
#include <LcLinAlgebra.h>

namespace Lucee
{
// constants to define integration method to use
  const static unsigned INT_RK4 = 0;
  const static unsigned INT_SEMI_IMPLICIT = 1;

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

// check if integration method specified
    intMethod = INT_RK4;
    if (tbl.hasString("method"))
    {
      std::string method = tbl.getString("method");
      if (method == "rk4")
        intMethod = INT_RK4;
      else if (method == "semi-implicit")
        intMethod = INT_SEMI_IMPLICIT;
      else
      {
        Lucee::Except lde("GridOdePointIntegrator::readInput: Unknown integration method '");
        lde << method << "' specified." << std::endl;
        throw lde;
      }
    }
  }

  template <unsigned NDIM>
  void
  GridOdePointIntegrator<NDIM>::integrate(double t0, double t1, Lucee::Field<NDIM, double>& sol)
  {
    ts.resize(sol.getNumComponents());
    switch (intMethod)
    {
      case INT_RK4:
        rk4(t0, t1-t0, sol);
        break;

      case INT_SEMI_IMPLICIT:
        semiImplicit(t0, t1-t0, sol);
        break;

      default:
// can not happen as all possible cases are tested in readInput() method
        break;
    }
  }

  template <unsigned NDIM>
  void
  GridOdePointIntegrator<NDIM>::rk4(double t0, double dt, Lucee::Field<NDIM, double>& sol)
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
      seq.fillWithIndex(idx);
// get centroid coordinate
      grid.setIndex(idx);
      grid.getCentroid(xc);
// set pointers
      sol.setPtr(inpPtr, idx);

// RK stage 1
      calcSource(t0, xc, &inpPtr[0], src);
      for (unsigned i=0; i<n; ++i)
        ql[i] = inpPtr[i] + hh*src[i];

// RK stage 2
      calcSource(t0+hh, xc, &ql[0], srct);
      for (unsigned i=0; i<n; ++i)
        ql[i] = inpPtr[i] + hh*srct[i];

// RK stage 3
      calcSource(t0+hh, xc, &ql[0], srcm);
      for (unsigned i=0; i<n; ++i)
      {
        ql[i] = inpPtr[i] + dt*srcm[i];
        srcm[i] = srct[i] + srcm[i];
      }

// RK stage 4
      calcSource(t0+dt, xc, &ql[0], srct);
// perform final update
      sol.setPtr(solPtr, idx);
      for (unsigned i=0; i<n; ++i)
        solPtr[i] += h6*(src[i] + srct[i] + 2*srcm[i]);
    }
  }

  template <unsigned NDIM>
  void
  GridOdePointIntegrator<NDIM>::semiImplicit(double t0, double dt, Lucee::Field<NDIM, double>& sol)
  {
// number of components to update
    unsigned n = sol.getNumComponents();
// allocate memory for storing source and Jacobian
    std::vector<double> src(n);
    Lucee::Matrix<double> unit(n, n), srcJac(n, n), rhs(n, 1);

// initialize unit matrix
    unit = 0.0;
    for (unsigned i=n; i<n; ++i)
      unit(i,i) = 1.0;

// get reference to grid
    const Lucee::StructuredGridBase<NDIM>& grid = this->getGrid();

// create pointers to fields
    Lucee::FieldPtr<double> solPtr = sol.createPtr();
    Lucee::ConstFieldPtr<double> inpPtr = sol.createConstPtr();

    double hh = dt/2.0;

    int idx[NDIM];
    double xc[3];
    Lucee::RowMajorSequencer<NDIM> seq(sol.getRegion());
// loop over field, updating solution in each cell
    while (seq.step())
    {
      seq.fillWithIndex(idx);
// get centroid coordinate
      grid.setIndex(idx);
      grid.getCentroid(xc);
// set pointers
      sol.setPtr(inpPtr, idx);

// compute compute source terms
      calcSource(t0, xc, &inpPtr[0], src);
// compute source Jacobian
      calcSourceJac(t0, xc, &inpPtr[0], srcJac);
// now compute matrix to invert
      for (unsigned i=0; i<n; ++i)
      {
        rhs(i,0) = src[i];
        for (unsigned j=0; j<n; ++j)
          srcJac(i,j) = unit(i,j)-hh*srcJac(i,j);
      }
// invert matrix and multiply it by source
      Lucee::solve(srcJac, rhs);
// perform final update
      sol.setPtr(solPtr, idx);
      for (unsigned i=0; i<n; ++i)
        solPtr[i] += dt*rhs(i,1);
    }

  }

  template <unsigned NDIM>
  void
  GridOdePointIntegrator<NDIM>::calcSource(double tm, const double xc[3], const double *inp, std::vector<double>& src)
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
      rhs[i]->calcSource(tm, xc, inp, &ts[0]);
      for (unsigned k=0; k<n; ++k)
        src[k] += ts[k];
    }
  }

  template <unsigned NDIM>
  void
  GridOdePointIntegrator<NDIM>::calcSourceJac(double tm, const double xc[3], 
    const double *inp, Lucee::Matrix<double>& srcJac)
  {
    unsigned n = srcJac.numRows();
    Lucee::Matrix<double> myJac(n, n);
// zap source jacobian
    srcJac = 0.0;
// compute each jacobian, accumulating it
    for (unsigned i=0; i<rhs.size(); ++i)
    {
      myJac = 0.0;
      rhs[i]->calcSourceJac(tm, xc, inp, myJac);
      for (unsigned ix=0; ix<n; ++ix)
        for (unsigned iy=0; iy<n; ++iy)
          srcJac(ix,iy) += myJac(ix,iy);
    }
  }

// instantiations
  template class GridOdePointIntegrator<1>;
  template class GridOdePointIntegrator<2>;
  template class GridOdePointIntegrator<3>;
}
