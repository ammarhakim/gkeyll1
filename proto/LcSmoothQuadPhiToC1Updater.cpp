/**
 * @file	LcSmoothQuadPhiToC1Updater.cpp
 *
 * @brief	Project piecewise quadriatic phi to C1 basis functions
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLinAlgebra.h>
#include <LcSmoothQuadPhiToC1Updater.h>

namespace Lucee
{
  static const char *id;

  SmoothQuadPhiToC1Updater::SmoothQuadPhiToC1Updater()
    : UpdaterIfc()
  {
  }

  void
  SmoothQuadPhiToC1Updater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis");
    else
      throw Lucee::Except("SmoothQuadPhiToC1Updater::readInput: Must specify element to use using 'basis'");
  }

  void
  SmoothQuadPhiToC1Updater::initialize()
  {
    Lucee::UpdaterIfc::initialize();

// get hold of grid
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();

// local region to update
    Lucee::Region<1, int> localRgn = grid.getLocalRegion();
    Lucee::RowMajorSequencer<1> seq(localRgn);
    seq.step();
    int idx[1];
    seq.fillWithIndex(idx);

// set index to first location in grid (this is okay as in this
// updater we are assuming grid is uniform)
    nodalBasis->setIndex(idx);

    unsigned nVol = nodalBasis->getNumGaussNodes();
    unsigned nlocal = nodalBasis->getNumNodes();

// get stiffness matrix
    Lucee::Matrix<double> stiffMatrix(nlocal, nlocal);
    nodalBasis->getGradStiffnessMatrix(0, stiffMatrix);

// calculate differentiation matrix
    diffMatrix = Lucee::Matrix<double>(nlocal, nlocal);
    for (unsigned i=0; i<nlocal; ++i)
      for (unsigned j=0; j<nlocal; ++j)
// diff matrices are computed from transposed stiffness matrices
          diffMatrix(i,j) = stiffMatrix(j,i);

// multiply by inverse of mass matrix
    Lucee::Matrix<double> massMatrix(nlocal, nlocal);
    nodalBasis->getMassMatrix(massMatrix);
    Lucee::solve(massMatrix, diffMatrix);
  }

  Lucee::UpdaterStatus
  SmoothQuadPhiToC1Updater::update(double t)
  {

    return Lucee::UpdaterStatus();
  }

  void
  SmoothQuadPhiToC1Updater::declareTypes()
  {
// takes one input (quadriatic, C0 phi)
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
// returns one output, C1 phi
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }
}
