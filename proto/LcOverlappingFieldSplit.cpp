/**
 * @file	LcOverlappingFieldSplit.cpp
 *
 * @brief	Updater to copy data from global domain to two split domains
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcOverlappingFieldSplit.h>
#include <LcStructuredGridBase.h>

namespace Lucee
{
  template <> const char *OverlappingFieldSplit<1>::id = "OverlappingFieldSplit1D";
  template <> const char *OverlappingFieldSplit<2>::id = "OverlappingFieldSplit2D";
  template <> const char *OverlappingFieldSplit<3>::id = "OverlappingFieldSplit3D";

  template <unsigned NDIM>
  OverlappingFieldSplit<NDIM>::OverlappingFieldSplit()
    : Lucee::UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  OverlappingFieldSplit<NDIM>::~OverlappingFieldSplit()
  {
  }

  template <unsigned NDIM>
  void
  OverlappingFieldSplit<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method first
    Lucee::UpdaterIfc::readInput(tbl);

    numCells = tbl.getNumber("numOverlappingCells");
    dir = tbl.getNumber("dir");
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  OverlappingFieldSplit<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    const Lucee::Field<NDIM, double>& qFull = this->getInp<Lucee::Field<NDIM, double> >(0);
    
    Lucee::Field<NDIM, double>& qLeft = this->getOut<Lucee::Field<NDIM, double> >(0);;
    Lucee::Field<NDIM, double>& qRight = this->getOut<Lucee::Field<NDIM, double> >(1);

    int idxF[NDIM], idxL[NDIM], idxR[NDIM];

    Lucee::ConstFieldPtr<double> ptrF = qFull.createConstPtr();
    Lucee::FieldPtr<double> ptrL = qLeft.createPtr();
    Lucee::FieldPtr<double> ptrR = qRight.createPtr();

    Lucee::Region<NDIM, int> leftRgn = qLeft.getRegion();
    Lucee::Region<NDIM, int> rightRgn = qRight.getRegion();
    Lucee::Region<NDIM, int> fullRgn = qFull.getRegion();

    int overlapLeftIdx = fullRgn.getUpper(dir) - rightRgn.getShape(dir);
    int overlapRightIdx = fullRgn.getLower(dir) + leftRgn.getShape(dir);

    Lucee::RowMajorSequencer<NDIM> seqFull(qFull.getExtRegion());
    while (seqFull.step())
    {
      seqFull.fillWithIndex(idxF);
      seqFull.fillWithIndex(idxL);
      seqFull.fillWithIndex(idxR);

      qFull.setPtr(ptrF, idxF);

      qLeft.setPtr(ptrL, idxL);
      idxR[dir] = idxF[dir]-overlapLeftIdx; // adjust index
      qRight.setPtr(ptrR, idxR);

      if (idxF[dir] < overlapLeftIdx)
      { // left region
        for (unsigned k=0; k<qFull.getNumComponents(); ++k)
          ptrL[k] = ptrF[k];

        // set left ghost cell of right field
        if (idxF[dir] == overlapLeftIdx-1)
        {
          for (unsigned k=0; k<qFull.getNumComponents(); ++k)
            ptrR[k] = ptrF[k];
        }
      }
      else if (idxF[dir] >= overlapRightIdx)
      { // right region
        for (unsigned k=0; k<qFull.getNumComponents(); ++k)
          ptrR[k] = ptrF[k];

        // set right ghost cell of left field
        if (idxF[dir] == overlapRightIdx)
        {
          for (unsigned k=0; k<qFull.getNumComponents(); ++k)
            ptrL[k] = ptrF[k];
        }        
      }
      else
      { // overlap region
        for (unsigned k=0; k<qFull.getNumComponents(); ++k)
        {
          ptrL[k] = ptrF[k];
          ptrR[k] = ptrF[k];
        }
      }
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  OverlappingFieldSplit<NDIM>::declareTypes()
  {
// expect one input field
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
// expect two output field
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));    
  }

// instantiations
  template class OverlappingFieldSplit<1>;
  template class OverlappingFieldSplit<2>;
  template class OverlappingFieldSplit<3>;
}
