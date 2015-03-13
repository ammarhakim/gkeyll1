/**
 * @file	LcEnstrophyUpdater.h
 *
 * @brief	Updater to compute energy from streamfunction.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcEnstrophyUpdater.h>
#include <LcGlobals.h>
#include <LcLinAlgebra.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <vector>

namespace Lucee
{
  template <> const char *EnstrophyUpdater<1>::id = "TotalEnstrophy1D";
  template <> const char *EnstrophyUpdater<2>::id = "TotalEnstrophy"; // NOTE: this is for backward compatibility
  template <> const char *EnstrophyUpdater<3>::id = "TotalEnstrophy3D";

  template <unsigned NDIM>
  EnstrophyUpdater<NDIM>::EnstrophyUpdater()
    : Lucee::UpdaterIfc()
  {
  }
  
  template <unsigned NDIM>
  void
  EnstrophyUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::UpdaterIfc::readInput(tbl);

// get hold of element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("NodalFiniteElementIfc::readInput: Must specify element to use using 'basis'");
  }

  template <unsigned NDIM>
  void
  EnstrophyUpdater<NDIM>::initialize()
  {
// call base class method
    Lucee::UpdaterIfc::initialize();

// get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

// local region to update
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    int idx[NDIM];
    for (unsigned d=0; d<NDIM; ++d)
      idx[d] = localRgn.getLower(d);
// set index to first location in grid (this is okay as in this
// updater we are assuming grid is uniform)
    nodalBasis->setIndex(idx);

    unsigned nVol = nodalBasis->getNumGaussNodes();
    unsigned nlocal = nodalBasis->getNumNodes();

// allocate space to get Gaussian quadrature data
    interpMat.m = Lucee::Matrix<double>(nVol, nlocal);
    ordinates.m = Lucee::Matrix<double>(nVol, 3);
    weights.resize(nVol);

    nodalBasis->getGaussQuadData(interpMat.m, ordinates.m, weights);
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  EnstrophyUpdater<NDIM>::update(double t)
  {
// get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
// get input arrays
    const Lucee::Field<NDIM, double>& chi = this->getInp<Lucee::Field<NDIM, double> >(0);
// get output dynVector
    Lucee::DynVector<double>& enstrophy = this->getOut<Lucee::DynVector<double> >(0);

// local region to update
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

// number of ordinates for volume integral
    unsigned nVol = nodalBasis->getNumGaussNodes();
// number of local nodes
    unsigned nlocal = nodalBasis->getNumNodes();

    Lucee::ConstFieldPtr<double> chiPtr = chi.createConstPtr();

// space for various quantities
    std::vector<double> quadChi(nVol);

    int idx[NDIM];
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);

    double totalEnstrophy = 0.0;
// compute total enstrophy
    while (seq.step())
    {
      seq.fillWithIndex(idx);

      nodalBasis->setIndex(idx);
      chi.setPtr(chiPtr, idx);

// compute vorticity at quadrature nodes
      matVec(1.0, interpMat.m, &chiPtr[0], 0.0, &quadChi[0]);
// compute contribution to enstrophy from this cell
      for (unsigned k=0; k<nVol; ++k)
        totalEnstrophy += weights[k]*quadChi[k]*quadChi[k];
    }

    double netTotalEnstrophy = totalEnstrophy;
// get hold of comm pointer to do all parallel messaging
    TxCommBase *comm = this->getComm();
// sum across all processors
    comm->allreduce(1, &totalEnstrophy, &netTotalEnstrophy, TX_SUM);

    std::vector<double> data(1);
    data[0] = 0.5*netTotalEnstrophy;
// push value into dynVector
    enstrophy.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  EnstrophyUpdater<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }

  template <unsigned NDIM>
  void 
  EnstrophyUpdater<NDIM>::matVec(double m, const Lucee::Matrix<double>& mat,
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
  template class EnstrophyUpdater<1>;
  template class EnstrophyUpdater<2>;
  template class EnstrophyUpdater<3>;
}
