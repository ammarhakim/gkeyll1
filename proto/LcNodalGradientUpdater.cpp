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
  template <> const char *NodalGradientUpdater<1>::id = "NodalGradient1D";
  template <> const char *NodalGradientUpdater<2>::id = "NodalGradient2D";
  template <> const char *NodalGradientUpdater<3>::id = "NodalGradient3D";

  template <unsigned NDIM>
  NodalGradientUpdater<NDIM>::NodalGradientUpdater()
    : Lucee::UpdaterIfc()
  {
  }
  
  template <unsigned NDIM>
  void
  NodalGradientUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::UpdaterIfc::readInput(tbl);

// get hold of element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("NodalGradientUpdater::readInput: Must specify element to use using 'basis'");
  }

  template <unsigned NDIM>
  void
  NodalGradientUpdater<NDIM>::initialize()
  {
// call base class method
    Lucee::UpdaterIfc::initialize();

// get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

// local region to update
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    seq.step();
    int idx[NDIM];
    seq.fillWithIndex(idx);

// set index to first location in grid (this is okay as in this
// updater we are assuming grid is uniform)
    nodalBasis->setIndex(idx);

    unsigned nVol = nodalBasis->getNumGaussNodes();
    unsigned nlocal = nodalBasis->getNumNodes();

// space for mass matrix
    Lucee::Matrix<double> massMatrix(nlocal, nlocal);

    for (unsigned dir=0; dir<NDIM; ++dir)
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

    for (unsigned dir=0; dir<NDIM; ++dir)
    {
      pDiffMatrix[dir].m = Lucee::Matrix<double>(nVol, nlocal);
// compute differentiation matrices that compute derivatives at
// quadrature nodes
      Lucee::accumulate(pDiffMatrix[dir].m, interpMat.m, diffMatrix[dir].m);
    }
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  NodalGradientUpdater<NDIM>::update(double t)
  {
// get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
// get input arrays
    const Lucee::Field<NDIM, double>& phi = this->getInp<Lucee::Field<NDIM, double> >(0);
// get output arrays
    Lucee::Field<NDIM, double>& gradXY = this->getOut<Lucee::Field<NDIM, double> >(0);

// local region to update
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);

// number of nodes in quadrature
    unsigned nVol = nodalBasis->getNumGaussNodes();
// number of local nodes
    unsigned nlocal = nodalBasis->getNumNodes();

// space for various quantities
    std::vector<double> grad(nlocal);
    Lucee::FieldPtr<double> gPtr = gradXY.createPtr();
    Lucee::ConstFieldPtr<double> phiK = phi.createConstPtr();

    int idx[NDIM];
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      nodalBasis->setIndex(idx);
      gradXY.setPtr(gPtr, idx);
      phi.setPtr(phiK, idx);

// extract potential at this location
      //nodalBasis->extractFromField(phi, phiK);

      for (unsigned dir=0; dir<NDIM; ++dir)
      {
// compute gradients in dir directions
        matVec(1.0, diffMatrix[dir].m, &phiK[0], 0.0, &grad[0]);
// copy gradient into output field
        for (unsigned k=0; k<nlocal; ++k)
          gPtr[NDIM*k+dir] = grad[k];
      }
    }
    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  NodalGradientUpdater<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void 
  NodalGradientUpdater<NDIM>::matVec(double m, const Lucee::Matrix<double>& mat,
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

// instantiations
  template class Lucee::NodalGradientUpdater<1>;
  template class Lucee::NodalGradientUpdater<2>;
  template class Lucee::NodalGradientUpdater<3>;
}
