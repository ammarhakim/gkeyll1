/**
 * @file	LcCopyContToDisContFieldUpdater.cpp
 *
 * @brief	Updater to copy continuous field to a discontinuous field.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCopyContToDisContFieldUpdater.h>
#include <LcField.h>

namespace Lucee
{
  template <> const char *CopyContToDisContFieldUpdater<1>::id = "CopyContToDisCont1D";
  template <> const char *CopyContToDisContFieldUpdater<2>::id = "CopyContToDisCont2D";
  template <> const char *CopyContToDisContFieldUpdater<3>::id = "CopyContToDisCont3D";

  template <unsigned NDIM>
  CopyContToDisContFieldUpdater<NDIM>::CopyContToDisContFieldUpdater()
    : Lucee::UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  CopyContToDisContFieldUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except(
        "CopyContToDisContFieldUpdater::readInput: Must specify element to use using 'basis'");
  }

  template <unsigned NDIM>
  void
  CopyContToDisContFieldUpdater<NDIM>::initialize()
  {
// call base class method
    Lucee::UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  CopyContToDisContFieldUpdater<NDIM>::update(double t)
  {
// get input/output fields
    const Lucee::Field<NDIM, double>& cIn = this->getInp<Lucee::Field<NDIM, double> >(0);
    Lucee::Field<NDIM, double>& dcOut = this->getOut<Lucee::Field<NDIM, double> >(0); 

    unsigned nlocal = nodalBasis->getNumNodes();
    std::vector<double> cData(nlocal);

    int idx[NDIM];
    Lucee::FieldPtr<double> dcPtr = dcOut.createPtr();
// loop over each cell in extended region (but not ghost cells)
    Lucee::Region<NDIM, int> localExtRgn 
      = cIn.getExtRegion().intersect(cIn.getGlobalRegion());
    Lucee::RowMajorSequencer<NDIM> seq(localExtRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      nodalBasis->setIndex(idx);
      
// get data from continous field
      nodalBasis->extractFromField(cIn, cData);

// copy it into discontinuous field
      dcOut.setPtr(dcPtr, idx);
      for (unsigned k=0; k<nlocal; ++k)
        dcPtr[k] = cData[k];
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  CopyContToDisContFieldUpdater<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class CopyContToDisContFieldUpdater<1>;
  template class CopyContToDisContFieldUpdater<2>;
  template class CopyContToDisContFieldUpdater<3>;
}
