/**
 * @file	LcModalL2NormUpdater.cpp
 *
 * @brief	Updater to compute the l2 norm of a field.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcModalL2NormUpdater.h>
#include <LcGlobals.h>
#include <LcStructuredGridBase.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <vector>

namespace Lucee
{
  const char *ModalL2NormUpdater::id = "ModalL2Norm";

  ModalL2NormUpdater::ModalL2NormUpdater()
    : Lucee::UpdaterIfc()
  {
  }
  
  void
  ModalL2NormUpdater::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    // Get number of basis functions to project on
    numBasis = (unsigned) tbl.getNumber("numBasis");
  }

  void
  ModalL2NormUpdater::initialize()
  {
    // call base class method
    Lucee::UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus
  ModalL2NormUpdater::update(double t)
  {
    // get hold of grid
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();
    // get input/output arrays
    const Lucee::Field<1, double>& q = this->getInp<Lucee::Field<1, double> >(0);
    Lucee::DynVector<double>& l2value = this->getOut<Lucee::DynVector<double> >(0);

    // local region to index
    Lucee::Region<1, int> localRgn = grid.getLocalRegion();
    // iterators
    Lucee::ConstFieldPtr<double> qPtr = q.createConstPtr();

    double dx = grid.getDx(0);

    double totalField = 0.0;
    // compute total enstrophy
    for (int i = localRgn.getLower(0); i < localRgn.getUpper(0); i++)
    {
      q.setPtr(qPtr, i); // right cell
      // Compute f^2 in the cell
      for (int basisIndex = 0; basisIndex < numBasis; basisIndex++)
      {
        totalField += qPtr[basisIndex]*qPtr[basisIndex]*dx/(2*basisIndex+1);
      }
    }

    double netTotalField = totalField;
// get hold of comm pointer to do all parallel messaging
    TxCommBase *comm = this->getComm();
// sum across all processors
    comm->allreduce(1, &totalField, &netTotalField, TX_SUM);

    std::vector<double> data(1);
    data[0] = netTotalField;
// push value into dynVector
    l2value.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  void
  ModalL2NormUpdater::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }
}
