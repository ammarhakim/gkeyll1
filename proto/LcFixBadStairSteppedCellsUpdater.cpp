/**
 * @file	LcFixBadStairSteppedCellsUpdater.cpp
 *
 * @brief	Detect and fix degenerate cells in a stair-stepped mesh.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAlignedRectCoordSys.h>
#include <LcFixBadStairSteppedCellsUpdater.h>
#include <LcField.h>
#include <LcStructuredGridBase.h>

namespace Lucee
{
  template <> const char *FixBadStairSteppedCellsUpdater<2>::id = "FixBadStairSteppedCells2D";
  template <> const char *FixBadStairSteppedCellsUpdater<3>::id = "FixBadStairSteppedCells3D";

  template <unsigned NDIM>
  FixBadStairSteppedCellsUpdater<NDIM>::FixBadStairSteppedCellsUpdater() 
    : Lucee::UpdaterIfc() 
  {
  }

  template <unsigned NDIM>
  FixBadStairSteppedCellsUpdater<NDIM>::~FixBadStairSteppedCellsUpdater() 
  {
  }

  template <unsigned NDIM>
  void
  FixBadStairSteppedCellsUpdater<NDIM>::readInput(Lucee::LuaTable& tbl) 
  {
    UpdaterIfc::readInput(tbl);
  }

  template <unsigned NDIM>
  void
  FixBadStairSteppedCellsUpdater<NDIM>::initialize() 
  {
// call base class method
    UpdaterIfc::initialize();

  }
  
  template <unsigned NDIM>
  Lucee::UpdaterStatus
  FixBadStairSteppedCellsUpdater<NDIM>::update(double t) 
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

// in/out field to fix    
    Lucee::Field<NDIM, double>& inOut =  this->getOut<Lucee::Field<NDIM, double> >(0);
    Lucee::FieldPtr<double> p = inOut.createPtr();
    Lucee::FieldPtr<double> pl = inOut.createPtr();
    Lucee::FieldPtr<double> pr = inOut.createPtr();

    numBadCells = 0; // keep count of bad cells
    for (unsigned dir=0; dir<NDIM; ++dir)
    {
// create sequencer to loop over *each* 1D slice in 'dir' direction
      Lucee::RowMajorSequencer<NDIM> seq(localRgn.deflate(dir));

// lower and upper bounds of 1D slice
      int sliceLower = localRgn.getLower(dir);
      int sliceUpper = localRgn.getUpper(dir);

      while (seq.step())
      {
        int idx[NDIM], idxr[NDIM], idxl[NDIM];
        seq.fillWithIndex(idx);
        seq.fillWithIndex(idxr);
        seq.fillWithIndex(idxl);
// loop over each slice
        for (int i=sliceLower; i<sliceUpper; ++i)
        {
          idx[dir] = i;
          idxr[dir] = i+1;
          idxl[dir] = i-1;

// set pointers
          inOut.setPtr(p, idx);
          inOut.setPtr(pl, idxl);
          inOut.setPtr(pr, idxr);

// Check or bad cell: a cell is bad if its two neighbors are of the
// opposite sign.
          if (p[0]>0)
          {
            if ((pl[0]<0) && (pr[0]<0))
            {
              // bad cell detected
              numBadCells += 1;
            }
          }
          else if (p[0]<0)
          {
            if ((pl[0]>0) && (pr[0]>0))
            {
              // bad cell detected
              numBadCells += 1;
            }            
          }
        }
      }
    }
    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  FixBadStairSteppedCellsUpdater<NDIM>::declareTypes() 
  {
// a single output field
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void
  FixBadStairSteppedCellsUpdater<NDIM>::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
// call base class Lua methods
    UpdaterIfc::appendLuaCallableMethods(lfm);
    lfm.appendFunc("numBadCells", luaGetNumBadCells);
  }

  template <unsigned NDIM>  
  int
  FixBadStairSteppedCellsUpdater<NDIM>::luaGetNumBadCells(lua_State *L)
  {
    FixBadStairSteppedCellsUpdater *p
      = Lucee::PointerHolder<FixBadStairSteppedCellsUpdater>::getObj(L);
    lua_pushnumber(L, p->numBadCells);
    return 1;
  }

// instantiations
  template class FixBadStairSteppedCellsUpdater<2>;
  template class FixBadStairSteppedCellsUpdater<3>;
}
