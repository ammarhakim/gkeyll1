/**
 * @file	LcRunningAverageOfFieldCalc.cpp
 *
 * @brief	Using the past N inputs, computes average of field at each node
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGlobals.h>
#include <LcRunningAverageOfFieldCalc.h>
#include <LcStructuredGridBase.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  template <> const char *RunningAverageOfFieldCalc<3>::id = "RunningAverageOfFieldCalc3D";

  template <unsigned NDIM>
  RunningAverageOfFieldCalc<NDIM>::RunningAverageOfFieldCalc()
    : Lucee::UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  RunningAverageOfFieldCalc<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasNumber("sampleSize"))
      sampleSize = tbl.getNumber("sampleSize");
    else
      Lucee::Except lce("RunningAverageOfFieldCalc::readInput: Must provide a sampleSize value");
  }

  template <unsigned NDIM>
  void
  RunningAverageOfFieldCalc<NDIM>::initialize()
  {
    // call base class method
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    // get local region
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    int numLocalCells = localRgn.getVolume();

    localRecord.resize(numLocalCells);
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  RunningAverageOfFieldCalc<NDIM>::update(double t)
  {
    // get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    // get input field
    const Lucee::Field<NDIM, double>& fld = this->getInp<Lucee::Field<NDIM, double> >(0);
    // get output field
    Lucee::Field<NDIM, double>& fldAvg = this->getOut<Lucee::Field<NDIM, double> >(0);
    
    // local region
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    // index array for looping over region
    int idx[NDIM];
    // pointers to input and output fields
    Lucee::ConstFieldPtr<double> fldPtr = fld.createConstPtr();
    Lucee::FieldPtr<double> fldAvgPtr = fldAvg.createPtr();
    
    Lucee::RowMajorSequencer<NDIM> localSeq = RowMajorSequencer<NDIM>(localRgn);
    Lucee::RowMajorIndexer<NDIM> localIdxr = RowMajorIndexer<NDIM>(localRgn);

    int nlocal = fld.getNumComponents();

    // clear out output field
    fldAvg = 0.0;

    while (localSeq.step())
    {
      localSeq.fillWithIndex(idx);
      int cellIndex = localIdxr.getIndex(idx);
      
      fld.setPtr(fldPtr, idx);
      fldAvg.setPtr(fldAvgPtr, idx);

      std::vector<double> cellValues(nlocal);

      for (int i = 0; i < nlocal; i++)
        cellValues[i] = fldPtr[i];

      // Add the new element of fldAvg to the localRecord structure
      localRecord[cellIndex].push_back(cellValues);

      // If record has more entires than maximum allowed, delete the oldest entry
      if (localRecord[cellIndex].size() > sampleSize)
        localRecord[cellIndex].pop_front();

      // Iterate over list to perform a time average
      for (std::list<std::vector<double> >::const_iterator it = localRecord[cellIndex].begin(),
          end = localRecord[cellIndex].end();
          it != end; it++)
      {
        // At this particular time sample, sum all nodes in this cell to output field
        for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
          fldAvgPtr[nodeIndex] += (*it)[nodeIndex];
      }

      // Divide fldAvgPtr by number of samples used in sum
      for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
        fldAvgPtr[nodeIndex] = fldAvgPtr[nodeIndex]/localRecord[cellIndex].size();
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  RunningAverageOfFieldCalc<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class Lucee::RunningAverageOfFieldCalc<3>;
}
