/**
 * @file	LcNodalVlasovUpdater.cpp
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
#include <LcNodalVlasovUpdater.h>
#include <LcStructuredGridBase.h>

// std includes
#include <vector>

namespace Lucee
{
// set id for module system
  template <> const char *NodalVlasovUpdater<1,1>::id = "NodalVlasov1X1V";
  template <> const char *NodalVlasovUpdater<1,2>::id = "NodalVlasov1X2V";
  template <> const char *NodalVlasovUpdater<1,3>::id = "NodalVlasov1X3V";
  template <> const char *NodalVlasovUpdater<2,2>::id = "NodalVlasov2X2V";
  template <> const char *NodalVlasovUpdater<2,3>::id = "NodalVlasov2X3V";
  template <> const char *NodalVlasovUpdater<3,3>::id = "NodalVlasov3X3V";

  template <unsigned CDIM, unsigned VDIM>
  NodalVlasovUpdater<CDIM,VDIM>::NodalVlasovUpdater()
    : UpdaterIfc()
  {
  }

  template <unsigned CDIM, unsigned VDIM>  
  void 
  NodalVlasovUpdater<CDIM,VDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<CDIM+VDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<CDIM+VDIM> >("basis");
    else
      throw Lucee::Except("NodalVlasovUpdater::readInput: Must specify element to use using 'basis'");

// directions to update
    if (tbl.hasNumVec("updateDirections"))
    {
      std::vector<double> ud = tbl.getNumVec("updateDirections");
      for (unsigned i=0; i<std::min<unsigned>(3, ud.size()); ++i)
      {
        unsigned d = (unsigned) ud[i];
        if (d<3)
          updateDims.push_back(d);
        else
          throw Lucee::Except("updateDirections must be a table with 0, 1, or 2");
      }
    }
    else
    {
      for (unsigned i=0; i<CDIM+VDIM; ++i)
        updateDims.push_back(i);
    }

    cfl = tbl.getNumber("cfl");
    cflm = 1.1*cfl; // use slightly large max CFL to avoid thrashing around

    onlyIncrement = false;
// when onlyIncrement flag is set contribution is not added to the
// input field, i.e. only increment is computed
    if (tbl.hasBool("onlyIncrement"))
      onlyIncrement = tbl.getBool("onlyIncrement");
  }

  template <unsigned CDIM, unsigned VDIM>
  void 
  NodalVlasovUpdater<CDIM,VDIM>::initialize()
  {
    Lucee::UpdaterIfc::initialize();
  }

  template <unsigned CDIM, unsigned VDIM>
  Lucee::UpdaterStatus 
  NodalVlasovUpdater<CDIM,VDIM>::update(double t)
  {
    return Lucee::UpdaterStatus();
  }

  template <unsigned CDIM, unsigned VDIM>  
  void
  NodalVlasovUpdater<CDIM,VDIM>::declareTypes()
  {
// distribution function
    this->appendInpVarType(typeid(Lucee::Field<CDIM+VDIM, double>));
// E and B field in a single field
    this->appendInpVarType(typeid(Lucee::Field<CDIM, double>));
// returns one output
    this->appendOutVarType(typeid(Lucee::Field<CDIM+VDIM, double>));
  }

  template <unsigned CDIM, unsigned VDIM>
  void 
  NodalVlasovUpdater<CDIM,VDIM>::matVec(double mc, const Lucee::Matrix<double>& mat,
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
  template class NodalVlasovUpdater<1,1>;
  template class NodalVlasovUpdater<1,2>;
  template class NodalVlasovUpdater<1,3>;
  template class NodalVlasovUpdater<2,2>;
  template class NodalVlasovUpdater<2,3>;
  template class NodalVlasovUpdater<3,3>;
}
