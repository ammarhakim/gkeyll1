/**
 * @file	LcCopyNodalFieldsUpdater.cpp
 *
 * @brief	Updater to copy one nodal field to another
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAlignedRectCoordSys.h>
#include <LcField.h>
#include <LcLinAlgebra.h>
#include <LcMathLib.h>
#include <LcCopyNodalFields.h>
#include <LcStructuredGridBase.h>

// std includes
#include <vector>

namespace Lucee
{
// set id for module system
  template <> const char *CopyNodalFieldsUpdater<1,2>::id = "CopyNodalFields1D_2D";
  template <> const char *CopyNodalFieldsUpdater<1,3>::id = "CopyNodalFields1D_3V";
  template <> const char *CopyNodalFieldsUpdater<1,4>::id = "CopyNodalFields1D_4V";
  template <> const char *CopyNodalFieldsUpdater<2,4>::id = "CopyNodalFields2D_4D";
  template <> const char *CopyNodalFieldsUpdater<2,5>::id = "CopyNodalFields2D_5D";

  template <unsigned CDIM, unsigned PDIM>
  bool
  CopyNodalFieldsUpdater<CDIM,PDIM>::sameConfigCoords(unsigned n, unsigned cn, double dxMin,
    const Lucee::Matrix<double>& phaseC, const Lucee::Matrix<double>& confC)
  {
    for (unsigned d=0; d<CDIM; ++d)
      if (! (std::fabs(phaseC(n,d)-confC(cn,d))<1e-4*dxMin) )
        return false;
    return true;
  }

  template <unsigned CDIM, unsigned PDIM>
  CopyNodalFieldsUpdater<CDIM,PDIM>::CopyNodalFieldsUpdater()
    : UpdaterIfc()
  {
  }

  template <unsigned CDIM, unsigned PDIM>  
  void 
  CopyNodalFieldsUpdater<CDIM,PDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<PDIM> >("targetBasis"))
      phaseBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<PDIM> >("targetBasis");
    else
      throw Lucee::Except("CopyNodalFieldsUpdater::readInput: Must specify target-basis using 'targetBasis'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<CDIM> >("sourceBasis"))
      confBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<CDIM> >("sourceBasis");
    else
      throw Lucee::Except("CopyNodalFieldsUpdater::readInput: Must specify source-space basis using 'sourceBasis'");
  }

  template <unsigned CDIM, unsigned PDIM>
  void 
  CopyNodalFieldsUpdater<CDIM,PDIM>::initialize()
  {
    Lucee::UpdaterIfc::initialize();

// get hold of grid
    const Lucee::StructuredGridBase<PDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<PDIM> >();
// local region to update
    Lucee::Region<PDIM, int> localRgn = grid.getLocalRegion();

    Lucee::RowMajorSequencer<PDIM> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[PDIM];
    seq.fillWithIndex(idx);
    phaseBasis->setIndex(idx);
    confBasis->setIndex(idx); // only first CDIM elements are used
    
    unsigned nlocal = phaseBasis->getNumNodes();

// compute mapping of phase-space nodes to configuration space
// nodes. The assumption here is that the node layout in phase-space
// and configuration space are such that each node in phase-space has
// exactly one node co-located with it in configuration space. No
// "orphan" phase-space node are allowed, and an exception is thrown
// if that occurs.
    phaseConfMap.resize(nlocal);
    Lucee::Matrix<double> phaseNodeCoords(phaseBasis->getNumNodes(), PNC);
    Lucee::Matrix<double> confNodeCoords(confBasis->getNumNodes(), CNC);

    double dxMin = grid.getDx(0);
    for (unsigned d=1; d<CDIM; ++d)
      dxMin = std::min(dxMin, grid.getDx(d));

    phaseBasis->getNodalCoordinates(phaseNodeCoords);
    confBasis->getNodalCoordinates(confNodeCoords);
    for (unsigned n=0; n<nlocal; ++n)
    {
      bool pcFound = false;
      for (unsigned cn=0; cn<nlocal; ++cn)
        if (sameConfigCoords(n, cn, dxMin, phaseNodeCoords, confNodeCoords))
        {
          phaseConfMap[n] = cn;
          pcFound = true;
          break;
        }
      if (!pcFound)
      {
        Lucee::Except lce(
          "CopyNodalFieldsUpdater::readInput: No matching source node for target node ");
        lce << n;
        throw lce;
      }
    }
  }

  template <unsigned CDIM, unsigned PDIM>
  Lucee::UpdaterStatus 
  CopyNodalFieldsUpdater<CDIM,PDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<PDIM>& grid
      = this->getGrid<Lucee::StructuredGridBase<PDIM> >();

    const Lucee::Field<CDIM, double>& q = this->getInp<Lucee::Field<CDIM, double> >(0);
    Lucee::Field<PDIM, double>& qNew = this->getOut<Lucee::Field<PDIM, double> >(0);

    return Lucee::UpdaterStatus();
  }

  template <unsigned CDIM, unsigned PDIM>  
  void
  CopyNodalFieldsUpdater<CDIM,PDIM>::declareTypes()
  {
    const unsigned NDIM = CDIM+PDIM;    
// distribution function
    this->appendInpVarType(typeid(Lucee::Field<CDIM, double>));
// returns one output: updated distribution function
    this->appendOutVarType(typeid(Lucee::Field<PDIM, double>));
  }

// instantiations
  template class CopyNodalFieldsUpdater<1,2>;
  template class CopyNodalFieldsUpdater<1,3>;
  template class CopyNodalFieldsUpdater<1,4>;
  template class CopyNodalFieldsUpdater<2,4>;
  template class CopyNodalFieldsUpdater<2,5>;
}
