/**
 * @file	LcNormGradPhiUpdater.cpp
 *
 * @brief	Updater to compute |grad.p|^2 integrated over the domain.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcNormGradPhiUpdater.h>
#include <LcGlobals.h>
#include <LcLinAlgebra.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <vector>

namespace Lucee
{
  template <> const char *NormGradPhiUpdater<1>::id = "NormGrad1D";
  template <> const char *NormGradPhiUpdater<2>::id = "NormGrad2D";
  template <> const char *NormGradPhiUpdater<3>::id = "NormGrad3D";

  template <unsigned NDIM>
  NormGradPhiUpdater<NDIM>::NormGradPhiUpdater()
    : Lucee::UpdaterIfc()
  {
  }
  
  template <unsigned NDIM>
  void
  NormGradPhiUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

// get hold of element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("NodalFiniteElementIfc::readInput: Must specify element to use using 'basis'");
  }

  template <unsigned NDIM>
  void
  NormGradPhiUpdater<NDIM>::initialize()
  {
    Lucee::UpdaterIfc::initialize();

    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

// set index to first location in grid (this is okay as in this
// updater we are assuming grid is uniform)
    int idx[NDIM];
    for (unsigned d=0; d<NDIM; ++d) 
      idx[d] = localRgn.getLower(d);
    nodalBasis->setIndex(idx);

    unsigned nlocal = nodalBasis->getNumNodes();
    Lucee::Matrix<double> massMatrix(nlocal, nlocal);

    for (unsigned dir=0; dir<NDIM; ++dir)
    {
// get stiffness matrix
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

    unsigned nVol = nodalBasis->getNumGaussNodes();

    interpMat.m = Lucee::Matrix<double>(nVol, nlocal);
    ordinates.m = Lucee::Matrix<double>(nVol, 3);
    weights.resize(nVol);

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
  NormGradPhiUpdater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    const Lucee::Field<NDIM, double>& phi = this->getInp<Lucee::Field<NDIM, double> >(0);
    Lucee::DynVector<double>& energy = this->getOut<Lucee::DynVector<double> >(0);

    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    unsigned nVol = nodalBasis->getNumGaussNodes();
    unsigned nlocal = nodalBasis->getNumNodes();

    std::vector<double> phiK(nlocal), normGradPhi(nVol);

    double totalEnergy = 0.0;
    int idx[NDIM];
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      nodalBasis->setIndex(idx);

      nodalBasis->extractFromField(phi, phiK);
      calcNormGrad(phiK, normGradPhi);

      for (unsigned qp=0; qp<nVol; ++qp)
        totalEnergy += weights[qp]*normGradPhi[qp];
    }

// sum across all processors
    double netTotalEnergy = totalEnergy;
    TxCommBase *comm = this->getComm();
    comm->allreduce(1, &totalEnergy, &netTotalEnergy, TX_SUM);

    std::vector<double> data(1);
    data[0] = netTotalEnergy;
    energy.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  NormGradPhiUpdater<NDIM>::calcNormGrad(std::vector<double>& phiK,
    std::vector<double>& normGradPhi)
  {
    for (unsigned i=0; i<normGradPhi.size(); ++i)
      normGradPhi[i] = 0.0;

    std::vector<double> grad(normGradPhi.size());
    for (unsigned d=0; d<NDIM; ++d)
    {
      matVec(1.0, pDiffMatrix[d].m, &phiK[0], 0.0, &grad[0]);
      for (unsigned i=0; i<normGradPhi.size(); ++i)
        normGradPhi[i] += grad[i]*grad[i];
    }
  }

  template <unsigned NDIM>
  void
  NormGradPhiUpdater<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }

  template <unsigned NDIM>
  void 
  NormGradPhiUpdater<NDIM>::matVec(double m, const Lucee::Matrix<double>& mat,
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
  template class NormGradPhiUpdater<1>;
  template class NormGradPhiUpdater<2>;
  template class NormGradPhiUpdater<3>;
}
