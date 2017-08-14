/**
 * @file	LcOverlappingFieldCopy.cpp
 *
 * @brief	Updater to copy data into ghost cells between two overlapping fields
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcOverlappingFieldCopy.h>
#include <LcStructuredGridBase.h>

namespace Lucee
{
  template <> const char *OverlappingFieldCopy<1>::id = "OverlappingFieldCopy1D";
  template <> const char *OverlappingFieldCopy<2>::id = "OverlappingFieldCopy2D";
  template <> const char *OverlappingFieldCopy<3>::id = "OverlappingFieldCopy3D";

  template <unsigned NDIM>
  OverlappingFieldCopy<NDIM>::OverlappingFieldCopy()
    : Lucee::UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  OverlappingFieldCopy<NDIM>::~OverlappingFieldCopy()
  {
  }

  template <unsigned NDIM>
  void
  OverlappingFieldCopy<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method first
    Lucee::UpdaterIfc::readInput(tbl);

    numOverlappingCells = tbl.getNumber("numOverlappingCells");
    dir = tbl.getNumber("dir");

    copyPeriodicDirs = false;
    if (tbl.hasBool("copyPeriodicDirs"))
      copyPeriodicDirs = tbl.getBool("copyPeriodicDirs");
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  OverlappingFieldCopy<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    Lucee::Field<NDIM, double>& qLeft = this->getOut<Lucee::Field<NDIM, double> >(0);;
    Lucee::Field<NDIM, double>& qRight = this->getOut<Lucee::Field<NDIM, double> >(1);

    int lo[NDIM], up[NDIM];
    int idxG[NDIM], idxI[NDIM];

    Lucee::FieldPtr<double> ptrL = qLeft.createPtr();
    Lucee::FieldPtr<double> ptrR = qRight.createPtr();
      
// 
// Step 1: copy stuff into "left" field ghost cells
//    
    
// create a region to represent ghost layer of "left" field
    for (unsigned i=0; i<NDIM; ++i)
    { // whole region, including extended region
      lo[i] = qLeft.getGlobalLowerExt(i);
      up[i] = qLeft.getGlobalUpperExt(i);
    }
// adjust region so it only indexes ghost cells
    lo[dir] = qLeft.getGlobalUpper(dir);
// region must be local to processor
    Lucee::Region<NDIM, int> gstRgn = qLeft.getExtRegion().intersect(
      Lucee::Region<NDIM, int>(lo, up));

// loop over ghost cells, copying stuff over
    Lucee::RowMajorSequencer<NDIM> seqGst(gstRgn);
    while (seqGst.step())
    {
      seqGst.fillWithIndex(idxG);
      seqGst.fillWithIndex(idxI);
      idxI[dir] = numOverlappingCells;

      qLeft.setPtr(ptrL, idxG);
      qRight.setPtr(ptrR, idxI);

      for (unsigned k=0; k<qLeft.getNumComponents(); ++k)
        ptrL[k] = ptrR[k];
    }

// 
// Step 2: copy stuff into "right" field ghost cells
//    
    
// create a region to represent ghost layer of "right" field
    for (unsigned i=0; i<NDIM; ++i)
    { // whole region, including extended region
      lo[i] = qRight.getGlobalLowerExt(i);
      up[i] = qRight.getGlobalUpperExt(i);
    }
// adjust region so it only indexes ghost cells
    up[dir] = qRight.getGlobalLower(dir);
// region must be local to processor
    gstRgn = qRight.getExtRegion().intersect(
      Lucee::Region<NDIM, int>(lo, up));

// loop over ghost cells, copying stuff over
    Lucee::RowMajorSequencer<NDIM> seqGst1(gstRgn);
    while (seqGst1.step())
    {
      seqGst1.fillWithIndex(idxG);
      seqGst1.fillWithIndex(idxI);
      idxI[dir] = qLeft.getUpper(dir)-numOverlappingCells-1;

      qRight.setPtr(ptrR, idxG);
      qLeft.setPtr(ptrL, idxI);

      for (unsigned k=0; k<qLeft.getNumComponents(); ++k)
        ptrR[k] = ptrL[k];
    }

// apply periodic BC is needed    
    if (copyPeriodicDirs)
    {
// 
// Step 1: copy stuff into "left" field ghost cells
//    
    
// create a region to represent ghost layer of "left" field
      for (unsigned i=0; i<NDIM; ++i)
      { // whole region, including extended region
        lo[i] = qLeft.getGlobalLowerExt(i);
        up[i] = qLeft.getGlobalUpperExt(i);
      }
// adjust region so it only indexes ghost cells
      up[dir] = qLeft.getGlobalLower(dir);
// region must be local to processor
      Lucee::Region<NDIM, int> gstRgn = qLeft.getExtRegion().intersect(
        Lucee::Region<NDIM, int>(lo, up));

// loop over ghost cells, copying stuff over
      Lucee::RowMajorSequencer<NDIM> seqGst(gstRgn);
      while (seqGst.step())
      {
        seqGst.fillWithIndex(idxG);
        seqGst.fillWithIndex(idxI);
        idxI[dir] = qRight.getUpper(dir)-1;

        qLeft.setPtr(ptrL, idxG);
        qRight.setPtr(ptrR, idxI);

        for (unsigned k=0; k<qLeft.getNumComponents(); ++k)
          ptrL[k] = ptrR[k];
      }

// 
// Step 2: copy stuff into "right" field ghost cells
//    
    
// create a region to represent ghost layer of "right" field
      for (unsigned i=0; i<NDIM; ++i)
      { // whole region, including extended region
        lo[i] = qRight.getGlobalLowerExt(i);
        up[i] = qRight.getGlobalUpperExt(i);
      }
// adjust region so it only indexes ghost cells
      lo[dir] = qRight.getGlobalUpper(dir);
// region must be local to processor
      gstRgn = qRight.getExtRegion().intersect(
        Lucee::Region<NDIM, int>(lo, up));

// loop over ghost cells, copying stuff over
      Lucee::RowMajorSequencer<NDIM> seqGst1(gstRgn);
      while (seqGst1.step())
      {
        seqGst1.fillWithIndex(idxG);
        seqGst1.fillWithIndex(idxI);
        idxI[dir] = 0;

        qRight.setPtr(ptrR, idxG);
        qLeft.setPtr(ptrL, idxI);

        for (unsigned k=0; k<qLeft.getNumComponents(); ++k)
          ptrR[k] = ptrL[k];
      }      
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  OverlappingFieldCopy<NDIM>::declareTypes()
  {
// expect two output fiels
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));    
  }

// instantiations
  template class OverlappingFieldCopy<1>;
  template class OverlappingFieldCopy<2>;
  template class OverlappingFieldCopy<3>;
}
