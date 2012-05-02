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
  const char *EnstrophyUpdater::id = "TotalEnstrophy";

  EnstrophyUpdater::EnstrophyUpdater()
    : Lucee::UpdaterIfc()
  {
  }
  
  void
  EnstrophyUpdater::readInput(Lucee::LuaTable& tbl)
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
  EnstrophyUpdater::initialize()
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

    nodalBasis->getGaussQuadData(interpMat.m, ordinates.m, weights);
  }

  Lucee::UpdaterStatus
  EnstrophyUpdater::update(double t)
  {
// get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();
// get input arrays
    const Lucee::Field<2, double>& chi = this->getInp<Lucee::Field<2, double> >(0);
// get output dynVector
    Lucee::DynVector<double>& enstrophy = this->getOut<Lucee::DynVector<double> >(0);

// local region to update
    Lucee::Region<2, int> localRgn = grid.getLocalRegion();

// number of local nodes
    unsigned nlocal = nodalBasis->getNumNodes();

    Lucee::ConstFieldPtr<double> chiPtr = chi.createConstPtr();

// space for various quantities
    std::vector<double> quadChi(nlocal);

    double totalEnstrophy = 0.0;
// compute total enstrophy
    for (int ix=localRgn.getLower(0); ix<localRgn.getUpper(0); ++ix)
    {
      for (int iy=localRgn.getLower(1); iy<localRgn.getUpper(1); ++iy)
      {
        nodalBasis->setIndex(ix, iy);
        chi.setPtr(chiPtr, ix, iy);

// compute vorticity at quadrature nodes
        matVec(1.0, interpMat.m, &chiPtr[0], 0.0, &quadChi[0]);

// compute contribition to enstrophy from this cell
        for (unsigned k=0; k<nlocal; ++k)
          totalEnstrophy += weights[k]*quadChi[k]*quadChi[k];
      }
    }

    double netTotalEnstrophy = totalEnstrophy;
// get hold of comm pointer to do all parallel messaging
    TxCommBase *comm = Loki::SingletonHolder<Lucee::Globals>
      ::Instance().comm;
// sum across all processors
    comm->allreduce(1, &totalEnstrophy, &netTotalEnstrophy, TX_SUM);

    std::vector<double> data(1);
    data[0] = 0.5*netTotalEnstrophy;
// push value into dynVector
    enstrophy.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  void
  EnstrophyUpdater::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }

  void 
  EnstrophyUpdater::matVec(double m, const Lucee::Matrix<double>& mat,
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
