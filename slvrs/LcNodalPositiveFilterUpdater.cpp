/**
 * @file	LcNodalPositiveFilterUpdater.cpp
 *
 * @brief	Updater to solve hyperbolic equations with nodal DG scheme.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAlignedRectCoordSys.h>
#include <LcField.h>
#include <LcLinAlgebra.h>
#include <LcMathLib.h>
#include <LcNodalPositiveFilterUpdater.h>
#include <LcStructuredGridBase.h>

// std includes
#include <vector>

namespace Lucee
{
// set id for module system
  template <> const char *NodalPositiveFilterUpdater<1>::id = "NodalPositiveFilter1D";
  template <> const char *NodalPositiveFilterUpdater<2>::id = "NodalPositiveFilter2D";
  template <> const char *NodalPositiveFilterUpdater<3>::id = "NodalPositiveFilter3D";

  template <unsigned NDIM>
  NodalPositiveFilterUpdater<NDIM>::NodalPositiveFilterUpdater()
    : UpdaterIfc()
  {
  }

  template <unsigned NDIM>  
  void 
  NodalPositiveFilterUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("NodalPositiveFilterUpdater::readInput: Must specify element to use using 'basis'");

    if (tbl.hasObject<Lucee::HyperEquation>("equation"))
      equation = &tbl.getObjectAsBase<Lucee::HyperEquation>("equation");
    else
    {
      Lucee::Except lce("NodalPositiveFilterUpdater::readInput: Must specify an equation to solve!");
      throw lce;
    }
  }

  template <unsigned NDIM>
  void 
  NodalPositiveFilterUpdater<NDIM>::initialize()
  {
    Lucee::UpdaterIfc::initialize();
// get weights for quadrature
    weights.resize(nodalBasis->getNumNodes());
    nodalBasis->getWeights(weights);
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus 
  NodalPositiveFilterUpdater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    const Lucee::Field<NDIM, double>& qNew = this->getInp<Lucee::Field<NDIM, double> >(0);
    Lucee::Field<NDIM, double>& qOld = this->getOut<Lucee::Field<NDIM, double> >(0);

    unsigned nlocal = nodalBasis->getNumNodes();
    unsigned meqn = equation->getNumEqns();

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>  
  void
  NodalPositiveFilterUpdater<NDIM>::declareTypes()
  {
// takes one input
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
// returns one output
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void 
  NodalPositiveFilterUpdater<NDIM>::matVec(double mc, const Lucee::Matrix<double>& mat,
    unsigned meqn, const double* vec, double v, double *out)
  {
    double tv;
    unsigned rows = mat.numRows(), cols = mat.numColumns();
    for (unsigned m=0; m<meqn; ++m)
    {
      for (unsigned i=0; i<rows; ++i)
      {
        tv = 0.0;
        for (unsigned j=0; j<cols; ++j)
          tv += mat(i,j)*vec[meqn*j+m];
        out[meqn*i+m] = mc*tv + v*out[meqn*i+m];
      }
    }
  }

// instantiations
  template class NodalPositiveFilterUpdater<1>;
  template class NodalPositiveFilterUpdater<2>;
  template class NodalPositiveFilterUpdater<3>;
}
