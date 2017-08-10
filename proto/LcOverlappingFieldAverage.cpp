/**
 * @file	LcOverlappingFieldAverage.cpp
 *
 * @brief	Updater to compute increment from source terms for use in DG schemes
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcOverlappingFieldAverage.h>
#include <LcStructuredGridBase.h>

namespace Lucee
{
  template <> const char *OverlappingFieldAverage<1>::id = "OverlappingFieldAverage1D";
  template <> const char *OverlappingFieldAverage<2>::id = "OverlappingFieldAverage2D";
  template <> const char *OverlappingFieldAverage<3>::id = "OverlappingFieldAverage3D";

  template <unsigned NDIM>
  OverlappingFieldAverage<NDIM>::OverlappingFieldAverage()
    : Lucee::UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  OverlappingFieldAverage<NDIM>::~OverlappingFieldAverage()
  {
  }

  template <unsigned NDIM>
  void
  OverlappingFieldAverage<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method first
    Lucee::UpdaterIfc::readInput(tbl);

    numCells = tbl.getNumber("numOverlappingCells");
    dir = tbl.getNumber("dir");
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  OverlappingFieldAverage<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    const Lucee::Field<NDIM, double>& qLeft = this->getInp<Lucee::Field<NDIM, double> >(0);
    const Lucee::Field<NDIM, double>& qRight = this->getInp<Lucee::Field<NDIM, double> >(1);

    Lucee::Field<NDIM, double>& qFull = this->getOut<Lucee::Field<NDIM, double> >(0);

    int idxF[NDIM], idxL[NDIM], idxR[NDIM];

    Lucee::ConstFieldPtr<double> ptrL = qLeft.createConstPtr();
    Lucee::ConstFieldPtr<double> ptrR = qRight.createConstPtr();
    Lucee::FieldPtr<double> ptrF = qFull.createPtr();

    Lucee::Region<NDIM, int> leftRgn = qLeft.getRegion();
    Lucee::Region<NDIM, int> rightRgn = qRight.getRegion();
    Lucee::Region<NDIM, int> fullRgn = qFull.getRegion();

    unsigned overlapLeftIdx = fullRgn.getUpper(dir) - rightRgn.getShape(dir);
    unsigned overlapRightIdx = fullRgn.getLower(dir) + leftRgn.getShape(dir);

    Lucee::RowMajorSequencer<NDIM> seqFull(qFull.getRegion());
    while (seqFull.step())
    {
      seqFull.fillWithIndex(idxF);
      seqFull.fillWithIndex(idxL);
      seqFull.fillWithIndex(idxR);

      qFull.setPtr(ptrF, idxF);

      if (idxF[dir] < overlapLeftIdx)
      { // left region
        qLeft.setPtr(ptrL, idxL);        

        for (unsigned k=0; k<qFull.getNumComponents(); ++k)
          ptrF[k] = ptrL[k];
      }
      else if (idxF[dir] >= overlapRightIdx)
      { // right region
        idxR[dir] = idxF[dir]-overlapLeftIdx; // adjust index
        qRight.setPtr(ptrR, idxR);

        for (unsigned k=0; k<qFull.getNumComponents(); ++k)
          ptrF[k] = ptrR[k];
      }
      else
      { // overlap region
        qLeft.setPtr(ptrL, idxL);
        idxR[dir] = idxF[dir]-overlapLeftIdx; // adjust index
        qRight.setPtr(ptrR, idxR);

        for (unsigned k=0; k<qFull.getNumComponents(); ++k)
          ptrF[k] = 0.5*(ptrL[k]+ptrR[k]); // average values
      }
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  OverlappingFieldAverage<NDIM>::declareTypes()
  {
// expect two input fields
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
// expect one output field
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class OverlappingFieldAverage<1>;
  template class OverlappingFieldAverage<2>;
  template class OverlappingFieldAverage<3>;
}
