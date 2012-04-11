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
    std::vector<double> weights(nlocal);

    double totalEnstrophy = 0.0;
// compute total enstrophy
    for (int ix=localRgn.getLower(0); ix<localRgn.getUpper(0); ++ix)
    {
      for (int iy=localRgn.getLower(1); iy<localRgn.getUpper(1); ++iy)
      {
        nodalBasis->setIndex(ix, iy);
        chi.setPtr(chiPtr, ix, iy);

// get quadrature weights
        nodalBasis->getWeights(weights);

// compute contribition to enstrophy from this cell
        for (unsigned k=0; k<nlocal; ++k)
          totalEnstrophy += weights[k]*chiPtr[k]*chiPtr[k];
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
}
