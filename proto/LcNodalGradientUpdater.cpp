/**
 * @file	LcNodalGradientUpdater.h
 *
 * @brief	Updater to compute energy from streamfunction.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcNodalGradientUpdater.h>
#include <LcGlobals.h>
#include <LcLinAlgebra.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <vector>

namespace Lucee
{
  const char *NodalGradientUpdater::id = "NodalGradient2D";

  NodalGradientUpdater::NodalGradientUpdater()
    : Lucee::UpdaterIfc()
  {
  }
  
  void
  NodalGradientUpdater::readInput(Lucee::LuaTable& tbl)
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
  NodalGradientUpdater::initialize()
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

    unsigned nVol = nodalBasis->getNumGaussNodes();
    unsigned nlocal = nodalBasis->getNumNodes();

// space for mass matrix
    Lucee::Matrix<double> massMatrix(nlocal, nlocal);

    for (unsigned dir=0; dir<2; ++dir)
    {
// get stiffness matrice
      Lucee::Matrix<double> stiffMatrix(nlocal, nlocal);
      nodalBasis->getGradStiffnessMatrix(dir, stiffMatrix);

// calculate differentiation matrix
      diffMatrix[dir].m = Lucee::Matrix<double>(nlocal, nlocal);
      for (unsigned i=0; i<nlocal; ++i)
        for (unsigned j=0; j<nlocal; ++j)
// diff matrices are computed from transposed stiffness matrices
          diffMatrix[dir].m(i,j) = stiffMatrix(j,i);

// multiply matrices by inverse of mass matrix
      nodalBasis->getMassMatrix(massMatrix);
      Lucee::solve(massMatrix, diffMatrix[dir].m);
    }

// allocate space to get Gaussian quadrature data
    interpMat.m = Lucee::Matrix<double>(nVol, nlocal);
    ordinates.m = Lucee::Matrix<double>(nVol, 3);
    weights.resize(nVol);

// get data needed for Gaussian quadrature
    nodalBasis->getGaussQuadData(interpMat.m, ordinates.m, weights);

    for (unsigned dir=0; dir<2; ++dir)
    {
      pDiffMatrix[dir].m = Lucee::Matrix<double>(nVol, nlocal);
// compute differentiation matrices that compute derivatives at
// quadrature nodes
      Lucee::accumulate(pDiffMatrix[dir].m, interpMat.m, diffMatrix[dir].m);
    }
  }

  Lucee::UpdaterStatus
  NodalGradientUpdater::update(double t)
  {
// get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();
// get input arrays
    const Lucee::Field<2, double>& phi = this->getInp<Lucee::Field<2, double> >(0);
// get output arrays
    Lucee::Field<2, double>& gradXY = this->getOut<Lucee::Field<2, double> >(0);

// local region to update
    Lucee::Region<2, int> localRgn = grid.getLocalRegion();

// number of nodes in quadrature
    unsigned nVol = nodalBasis->getNumGaussNodes();
// number of local nodes
    unsigned nlocal = nodalBasis->getNumNodes();

// space for various quantities
    std::vector<double> phiK(nlocal), gradX(nlocal), gradY(nlocal);
    Lucee::FieldPtr<double> gPtr = gradXY.createPtr();

// compute gradients
    for (int ix=localRgn.getLower(0); ix<localRgn.getUpper(0); ++ix)
    {
      for (int iy=localRgn.getLower(1); iy<localRgn.getUpper(1); ++iy)
      {
        nodalBasis->setIndex(ix, iy);
        gradXY.setPtr(gPtr, ix, iy);

// extract potential at this location
        nodalBasis->extractFromField(phi, phiK);
// compute gradients in X- and Y-directions
        matVec(1.0, pDiffMatrix[0].m, phiK, 0.0, &gradX[0]);
        matVec(1.0, pDiffMatrix[1].m, phiK, 0.0, &gradY[0]);

// copy gradient into output field
        for (unsigned k=0; k<nlocal; ++k)
        {
          gPtr[2*k+0] = gradX[k];
          gPtr[2*k+1] = gradY[k];
        }
      }
    }

    return Lucee::UpdaterStatus();
  }

  void
  NodalGradientUpdater::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }

  void 
  NodalGradientUpdater::matVec(double m, const Lucee::Matrix<double>& mat,
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
