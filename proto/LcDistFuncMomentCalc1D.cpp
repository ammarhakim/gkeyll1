/**
 * @file	LcDistFuncMomentCalc1D.h
 *
 * @brief	Updater to compute moments of distribution function.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDistFuncMomentCalc1D.h>
#include <LcGlobals.h>
#include <LcLinAlgebra.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <cmath>
#include <vector>

namespace Lucee
{
  const char *DistFuncMomentCalc1D::id = "DistFuncMomentCalc1D";

  DistFuncMomentCalc1D::DistFuncMomentCalc1D()
    : Lucee::UpdaterIfc()
  {
  }
  
  void
  DistFuncMomentCalc1D::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::UpdaterIfc::readInput(tbl);

// get hold of 2D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis2d"))
      nodalBasis2d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis2d");
    else
      throw Lucee::Except(
        "DistFuncMomentCalc1D::readInput: Must specify 2D element to use using 'basis2d'");

// get hold of 1D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis1d"))
      nodalBasis1d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis1d");
    else
      throw Lucee::Except(
        "DistFuncMomentCalc1D::readInput: Must specify 2D element to use using 'basis1d'");

// get moment to compute
    if (tbl.hasNumber("moment"))
    moment = (unsigned) tbl.getNumber("moment");
    else
      throw Lucee::Except(
        "DistFuncMomentCalc1D::readInput: Must specify moment using 'moment'");
  }

  void
  DistFuncMomentCalc1D::initialize()
  {
// call base class method
    Lucee::UpdaterIfc::initialize();

// get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();

// local region to update
    Lucee::Region<2, int> localRgn = grid.getLocalRegion();

// get number of nodes in 1D and 2D
    unsigned nlocal1d = nodalBasis1d->getNumNodes();
    unsigned nlocal2d = nodalBasis2d->getNumNodes();

// In the following the 2D element stores the moment matrix, while the
// 1D element stores the mass matrix (as it should). This is not
// really the correct way to do it: cross basis-function quantities
// should really be in their own class.

// allocate and get moment matrices
    for (unsigned m=0; m<3; ++m)
    {
      mm[m].m = Lucee::Matrix<double>(nlocal1d, nlocal2d);
      nodalBasis2d->getMomentMatrix(m, mm[m].m);
    }

    Lucee::Matrix<double> massMatrix1d(nlocal1d, nlocal1d);
// now multiply each of these by inverse of 1D mass matrix
    for (unsigned m=0; m<3; ++m)
    {
      nodalBasis1d->getMassMatrix(massMatrix1d);
      Lucee::solve(massMatrix1d, mm[m].m);
    }
  }

  Lucee::UpdaterStatus
  DistFuncMomentCalc1D::update(double t)
  {
// get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();

// get input field (2d)
    const Lucee::Field<2, double>& distF = this->getInp<Lucee::Field<2, double> >(0);
// get output field (1D)
    Lucee::Field<1, double>& moment = this->getOut<Lucee::Field<1, double> >(0);

// local region to update (This is the 2D region. The 1D region is
// assumed to have the same cell layout as the X-direction of the 2D region)
    Lucee::Region<2, int> localRgn = grid.getLocalRegion();

// clear out contents of output field
    moment = 0.0;

// iterators into fields
    Lucee::ConstFieldPtr<double> distFPtr = distF.createConstPtr();
    Lucee::FieldPtr<double> momentPtr = moment.createPtr();

// loop over all X-direction cells
    for (int i=localRgn.getLower(0); i<localRgn.getUpper(0); ++i)
    {
// set iterator into moment field (it is 1D)
      moment.setPtr(momentPtr, i);

// sum over all Y-direction cells
      for (int j=localRgn.getLower(1); j<localRgn.getUpper(1); ++j)
      {
// set iterator into distribution function (it is 2D)
        distF.setPtr(distFPtr, i, j);
// accumulate contribution to moment from this cell
        matVec(1.0, mm[0].m, &distFPtr[0], 1.0, &momentPtr[0]);
      }
    }

    return Lucee::UpdaterStatus();
  }

  void
  DistFuncMomentCalc1D::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }

  void 
  DistFuncMomentCalc1D::matVec(double m, const Lucee::Matrix<double>& mat,
    const double* vec, double v, double *out)
  {
    double tv;
    unsigned rows = mat.numRows(), cols = mat.numColumns();
    for (unsigned i=0; i<rows; ++i)
    {
      tv = 0.0;
      for (unsigned j=0; j<cols; ++j)
        tv += mat(i,j)*vec[j];
      out[i] = m*tv + v*out[i];
    }
  }
}
