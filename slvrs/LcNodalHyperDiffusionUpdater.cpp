/**
 * @file	LcNodalHyperDiffusionUpdater.cpp
 *
 * @brief	Updater to evaluate (hyper)diffusion operators using nodal DG
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLinAlgebra.h>
#include <LcNodalHyperDiffusionUpdater.h>


namespace Lucee
{
  template <> const char *NodalHyperDiffusionUpdater<1>::id = "HyperDiffusion1D";
  template <> const char *NodalHyperDiffusionUpdater<2>::id = "HyperDiffusion2D";
  template <> const char *NodalHyperDiffusionUpdater<3>::id = "HyperDiffusion3D";

  template <unsigned NDIM>
  NodalHyperDiffusionUpdater<NDIM>::NodalHyperDiffusionUpdater()
  {
  }

  template <unsigned NDIM>
  NodalHyperDiffusionUpdater<NDIM>::~NodalHyperDiffusionUpdater()
  {
  }

  template <unsigned NDIM>
  void
  NodalHyperDiffusionUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("NodalHyperDiffusionUpdater::readInput: Must specify element to use using 'basis'");

// diffusion coefficient
    alpha = tbl.getNumber("diffusionCoeff");
// CFL number to control time-step
    cfl = tbl.getNumber("cfl"); // CFL number
// should only increments be computed?
    onlyIncrement = false;
    if (tbl.hasBool("onlyIncrement"))
      onlyIncrement = tbl.getBool("onlyIncrement");

// directions to update
    if (tbl.hasNumVec("updateDirections"))
    {
      std::vector<double> ud = tbl.getNumVec("updateDirections");
      for (unsigned i=0; i<std::min<unsigned>(3, ud.size()); ++i)
      {
        unsigned d = (unsigned) ud[i];
        if (d<3)
          updateDirs.push_back(d);
        else
          throw Lucee::Except("updateDirections must be a table with 0, 1, or 2");
      }
    }
    else
    {
      for (unsigned i=0; i<NDIM; ++i)
        updateDirs.push_back(i);
    }
  }

  template <unsigned NDIM>
  void
  NodalHyperDiffusionUpdater<NDIM>::initialize()
  {
    Lucee::UpdaterIfc::initialize();
    
// get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
// local region to update
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    unsigned nlocal = nodalBasis->getNumNodes();
// allocate space for matrices
    iMat.resize(NDIM);
    lowerMat.resize(NDIM);
    upperMat.resize(NDIM);
    for (unsigned d=0; d<NDIM; ++d)
    {
      iMat[d] = Lucee::Matrix<double>(nlocal, nlocal);
      lowerMat[d] = Lucee::Matrix<double>(nlocal, nlocal);
      upperMat[d] = Lucee::Matrix<double>(nlocal, nlocal);
    }

    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[NDIM];
    seq.fillWithIndex(idx);
    nodalBasis->setIndex(idx);

// get various matrices needed
    nodalBasis->getDiffusionMatrices(iMat, lowerMat, upperMat);

// pre-multiply each of the matrices by inverse matrix
    Lucee::Matrix<double> massMatrix(nlocal, nlocal);

// NOTE: mass matrix is fetched repeatedly as the solve() method
// destroys it during the inversion process
    for (unsigned d=0; d<NDIM; ++d)
    {
      nodalBasis->getMassMatrix(massMatrix);
      Lucee::solve(massMatrix, iMat[d]);
    }
    for (unsigned d=0; d<NDIM; ++d)
    {
      nodalBasis->getMassMatrix(massMatrix);
      Lucee::solve(massMatrix, lowerMat[d]);
    }
    for (unsigned d=0; d<NDIM; ++d)
    {
      nodalBasis->getMassMatrix(massMatrix);
      Lucee::solve(massMatrix, upperMat[d]);
    }
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  NodalHyperDiffusionUpdater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    const Lucee::Field<NDIM, double>& inpFld = this->getInp<Lucee::Field<NDIM, double> >(0);
    Lucee::Field<NDIM, double>& diffOut = this->getOut<Lucee::Field<NDIM, double> >(0);

    double dt = t-this->getCurrTime();
    double dxMax = 0.0;
    for (unsigned d=0; d<NDIM; ++d)
      dxMax = std::max(dxMax, grid.getDx(d));

// check time-step
    double cflm = 1.1*cfl;
    double cfla = alpha*dt/(dxMax*dxMax);
    if (cfla>cflm)
      return Lucee::UpdaterStatus(false, dt*cfl/cfla);

    double fact = 1.0;
    if (onlyIncrement)
    {
// if only increments are requested, the updater computes alpha*d^2/dx^x inpFld
      diffOut = 0.0;
      fact = alpha;
    }
    else
    {
// updater computes inpFld + dt*alpha*d^2/dx^x inpFld
      diffOut.copy(inpFld);
      fact = alpha*dt;
    }
    
    Lucee::ConstFieldPtr<double> inpFldPtr = inpFld.createConstPtr();
    Lucee::FieldPtr<double> diffOutPtr = diffOut.createPtr();

    int idx[NDIM], idxL[NDIM];
// local region to index
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      diffOut.setPtr(diffOutPtr, idx);

      for (unsigned d=0; d<updateDirs.size(); ++d)
      {
        unsigned dir = updateDirs[d];
// add in contribution to cell from current cell
        inpFld.setPtr(inpFldPtr, idx);
        matVec(fact, iMat[d], &inpFldPtr[0], 1.0, &diffOutPtr[0]);

        seq.fillWithIndex(idxL);
// add in contribution from cells attached to lower/upper faces
        idxL[dir] = idx[dir]-1; // cell attached to lower face
        inpFld.setPtr(inpFldPtr, idxL);
        matVec(fact, lowerMat[d], &inpFldPtr[0], 1.0, &diffOutPtr[0]);

        idxL[dir] = idx[dir]+1; // cell attached to upper face
        inpFld.setPtr(inpFldPtr, idxL);
        matVec(fact, upperMat[d], &inpFldPtr[0], 1.0, &diffOutPtr[0]);
      }
    }

    return Lucee::UpdaterStatus(true, dt*cfl/cfla);
  }

  template <unsigned NDIM>
  void
  NodalHyperDiffusionUpdater<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void 
  NodalHyperDiffusionUpdater<NDIM>::matVec(double m, const Lucee::Matrix<double>& mat,
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
  template class NodalHyperDiffusionUpdater<1>;
  template class NodalHyperDiffusionUpdater<2>;
  template class NodalHyperDiffusionUpdater<3>;
}
