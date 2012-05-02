/**
 * @file	LcEnergyFromStreamFunctionUpdater.h
 *
 * @brief	Updater to compute energy from streamfunction.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcEnergyFromStreamFunctionUpdater.h>
#include <LcGlobals.h>
#include <LcLinAlgebra.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <vector>

namespace Lucee
{
  const char *EnergyFromStreamFunctionUpdater::id = "EnergyFromStreamFunction";

  EnergyFromStreamFunctionUpdater::EnergyFromStreamFunctionUpdater()
    : Lucee::UpdaterIfc()
  {
  }
  
  void
  EnergyFromStreamFunctionUpdater::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::UpdaterIfc::readInput(tbl);

// get hold of element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis");
    else
      throw Lucee::Except("NodalFiniteElementIfc::readInput: Must specify element to use using 'basis'");
  }

  void
  EnergyFromStreamFunctionUpdater::initialize()
  {
// call base class method
    Lucee::UpdaterIfc::initialize();

// get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();

// local region to update
    Lucee::Region<2, int> localRgn = grid.getLocalRegion();

// set index to first location in grid (this is okay as in this
// updater we are assuming grid is uniform)
    nodalBasis->setIndex(localRgn.getLower(0), localRgn.getLower(1));

    unsigned nlocal = nodalBasis->getNumNodes();

// allocate space to get Gaussian quadrature data
    interpMat.m = Lucee::Matrix<double>(nlocal, nlocal);
    ordinates.m = Lucee::Matrix<double>(nlocal, 3);
    weights.resize(nlocal);

    matVec(1.0, pDiffMatrix[0].m, phiK, 0.0, &gradX[0]);
    matVec(1.0, pDiffMatrix[1].m, phiK, 0.0, &gradY[0]);

// compute norm of gradient
    for (unsigned i=0; i<normGradPhi.size(); ++i)
      normGradPhi[i] = gradX[i]*gradX[i] + gradY[i]*gradY[i];
  }

  void
  EnergyFromStreamFunctionUpdater::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }

  void 
  EnergyFromStreamFunctionUpdater::matVec(double m, const Lucee::Matrix<double>& mat,
    const std::vector<double>& vec, double v, double *out)
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
