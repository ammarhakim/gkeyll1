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
  template <> const char *CopyNodalFieldsUpdater<2,3>::id = "CopyNodalFields2D_3D";
  template <> const char *CopyNodalFieldsUpdater<2,4>::id = "CopyNodalFields2D_4D";
  template <> const char *CopyNodalFieldsUpdater<2,5>::id = "CopyNodalFields2D_5D";
  template <> const char *CopyNodalFieldsUpdater<3,5>::id = "CopyNodalFields3D_5D";


  template <> const char *CopyNodalFieldsUpdater<1,1>::id = "CopyNodalFields1D"; 
  template <> const char *CopyNodalFieldsUpdater<2,2>::id = "CopyNodalFields2D";
  template <> const char *CopyNodalFieldsUpdater<3,3>::id = "CopyNodalFields3D";

  template <unsigned SDIM, unsigned TDIM>
  bool
  CopyNodalFieldsUpdater<SDIM,TDIM>::sameConfigCoords(unsigned n, unsigned cn, double dxMin,
    const Lucee::Matrix<double>& phaseC, const Lucee::Matrix<double>& confC)
  {
    for (unsigned d=0; d<SDIM; ++d)
      if (! (std::fabs(phaseC(n,coordinateMap[d])-confC(cn,d))<1e-4*dxMin) )
        return false;
    return true;
  }

  template <unsigned SDIM, unsigned TDIM>
  CopyNodalFieldsUpdater<SDIM,TDIM>::CopyNodalFieldsUpdater()
    : UpdaterIfc()
  {
  }

  template <unsigned SDIM, unsigned TDIM>  
  void 
  CopyNodalFieldsUpdater<SDIM,TDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<SDIM> >("sourceBasis"))
      sourceBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<SDIM> >("sourceBasis");
    else
      throw Lucee::Except("CopyNodalFieldsUpdater::readInput: Must specify source-space basis using 'sourceBasis'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<TDIM> >("targetBasis"))
      targetBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<TDIM> >("targetBasis");
    else
      throw Lucee::Except("CopyNodalFieldsUpdater::readInput: Must specify target-basis using 'targetBasis'");

    // Optional input for coordinate mapping. Otherwise, creates a vector
    // of size SDIM, representing the map (0,1,2,..) -> (0,1,2,..)
    if (tbl.hasNumVec("coordinateMap"))
    {
      std::vector<double> tempMap = tbl.getNumVec("coordinateMap");

      if(tempMap.size() != SDIM)
        throw Lucee::Except("CopyNodalFieldsUpdater::readInput: coordinateMap has an incorrect number of elements (!= SDIM)");

      for (int i = 0; i < SDIM; i++)
      {
        int d = (int) tempMap[i];
        if (d < TDIM)
          coordinateMap.push_back(d);
        else
          throw Lucee::Except("CopyNodalFieldsUpdater::readInput: coordinateMap must be a table with each element < TDIM");
      }
    }
    else
    {
      for (int i = 0; i < SDIM; i++)
        coordinateMap.push_back(i);
    }
  }

  template <unsigned SDIM, unsigned TDIM>
  void 
  CopyNodalFieldsUpdater<SDIM,TDIM>::initialize()
  {
    Lucee::UpdaterIfc::initialize();

// get hold of grid
    const Lucee::StructuredGridBase<TDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<TDIM> >();
// local region to update
    Lucee::Region<TDIM, int> localRgn = grid.getLocalRegion();

    Lucee::RowMajorSequencer<TDIM> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[TDIM];
    seq.fillWithIndex(idx);
    targetBasis->setIndex(idx);
    int idxSrc[SDIM];
    for (int sIndex = 0; sIndex < SDIM; sIndex++)
      idxSrc[sIndex] = idx[coordinateMap[sIndex]];
    sourceBasis->setIndex(idxSrc); // only first SDIM elements are used
    
    unsigned nlocal = targetBasis->getNumNodes();

// compute mapping of phase-space nodes to configuration space
// nodes. The assumption here is that the node layout in phase-space
// and configuration space are such that each node in phase-space has
// exactly one node co-located with it in configuration space. No
// "orphan" phase-space node are allowed, and an exception is thrown
// if that occurs.
    tarSrcMap.resize(nlocal);
    Lucee::Matrix<double> tarNodeCoords(targetBasis->getNumNodes(), PNC);
    Lucee::Matrix<double> srcNodeCoords(sourceBasis->getNumNodes(), CNC);

    double dxMin = grid.getDx(0);
    for (unsigned d=1; d<SDIM; ++d)
      dxMin = std::min(dxMin, grid.getDx(d));

    targetBasis->getNodalCoordinates(tarNodeCoords);
    sourceBasis->getNodalCoordinates(srcNodeCoords);

    for (unsigned n=0; n<nlocal; ++n)
    {
      bool pcFound = false;
      for (unsigned cn=0; cn<sourceBasis->getNumNodes(); ++cn)
        if (sameConfigCoords(n, cn, dxMin, tarNodeCoords, srcNodeCoords))
        {
          tarSrcMap[n] = cn;
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

  template <unsigned SDIM, unsigned TDIM>
  Lucee::UpdaterStatus 
  CopyNodalFieldsUpdater<SDIM,TDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<TDIM>& grid
      = this->getGrid<Lucee::StructuredGridBase<TDIM> >();

    const Lucee::Field<SDIM, double>& qSrc = this->getInp<Lucee::Field<SDIM, double> >(0);
    Lucee::Field<TDIM, double>& qTar = this->getOut<Lucee::Field<TDIM, double> >(0);

    Lucee::ConstFieldPtr<double> qSrcPtr = qSrc.createConstPtr();
    Lucee::FieldPtr<double> qTarPtr = qTar.createPtr();

    unsigned nlocal = targetBasis->getNumNodes();
    Lucee::Region<TDIM, int> localRgn = grid.getLocalRegion();
    int idx[TDIM];
    int idxSrc[SDIM];
    Lucee::RowMajorSequencer<TDIM> seq(localRgn);

    while (seq.step())
    {
      seq.fillWithIndex(idx);
      for (int sIndex = 0; sIndex < SDIM; sIndex++)
        idxSrc[sIndex] = idx[coordinateMap[sIndex]];
      qSrc.setPtr(qSrcPtr, idxSrc);
      qTar.setPtr(qTarPtr, idx);

      for (unsigned k=0; k<nlocal; ++k)
        qTarPtr[k] = qSrcPtr[tarSrcMap[k]];
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned SDIM, unsigned TDIM>  
  void
  CopyNodalFieldsUpdater<SDIM,TDIM>::declareTypes()
  {
// distribution function
    this->appendInpVarType(typeid(Lucee::Field<SDIM, double>));
// returns one output: updated distribution function
    this->appendOutVarType(typeid(Lucee::Field<TDIM, double>));
  }

// instantiations
  template class CopyNodalFieldsUpdater<1,2>;
  template class CopyNodalFieldsUpdater<1,3>;
  template class CopyNodalFieldsUpdater<1,4>;
  template class CopyNodalFieldsUpdater<2,3>;
  template class CopyNodalFieldsUpdater<2,4>;
  template class CopyNodalFieldsUpdater<2,5>;
  template class CopyNodalFieldsUpdater<3,5>;

  template class CopyNodalFieldsUpdater<1,1>;
  template class CopyNodalFieldsUpdater<2,2>;
  template class CopyNodalFieldsUpdater<3,3>;
}
