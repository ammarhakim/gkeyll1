/**
 * @file	LcStairSteppedBcUpdater.cpp
 *
 * @brief	Base class for boundary conditions for stair-stepped boundaries.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAlignedRectCoordSys.h>
#include <LcStairSteppedBcUpdater.h>
#include <LcField.h>
#include <LcStructuredGridBase.h>

namespace Lucee
{
// right/left flags
  enum {SSB_UP = 0, SSB_LO = 1};

  template <> const char *StairSteppedBcUpdater<1>::id = "StairSteppedBc1D";
  template <> const char *StairSteppedBcUpdater<2>::id = "StairSteppedBc2D";
  template <> const char *StairSteppedBcUpdater<3>::id = "StairSteppedBc3D";

  template <unsigned NDIM>
  StairSteppedBcUpdater<NDIM>::StairSteppedBcUpdater() 
    : Lucee::UpdaterIfc() 
  {
  }

  template <unsigned NDIM>
  StairSteppedBcUpdater<NDIM>::~StairSteppedBcUpdater() 
  {
    delete ssBnd;
  }

  template <unsigned NDIM>
  void
  StairSteppedBcUpdater<NDIM>::readInput(Lucee::LuaTable& tbl) 
  {
    UpdaterIfc::readInput(tbl);

    bcDir = 0;
// read direction to apply
    if (tbl.hasNumber("dir"))
      bcDir = tbl.getNumber("dir");

// get list of boundary conditions to apply
    if ( !tbl.hasTable("boundaryConditions") )
      Lucee::Except lce("StairSteppedBcUpdater::readInput: Must specify boundary conditions to apply!");

    Lucee::LuaTable bcTbl = tbl.getTable("boundaryConditions");
    bcList = bcTbl.template getAllObjects<Lucee::BoundaryCondition>();

// set pointer to in/out field
    inOut = &tbl.getObject<Lucee::Field<NDIM, double> >("inOutField");
  }

  template <unsigned NDIM>
  void
  StairSteppedBcUpdater<NDIM>::initialize() 
  {
// call base class method
    UpdaterIfc::initialize();

    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    int lg[NDIM], ug[NDIM];
    for (unsigned d=0; d<NDIM; ++d)
      lg[d] = ug[d] = 1; // just a single layer of ghost cells is needed
    ssBnd = new Lucee::Field<NDIM, double>(localRgn, NDIM, lg, ug); // to store boundary information

    Lucee::FieldPtr<double> iop = inOut->createPtr();
    Lucee::FieldPtr<double> iopr = inOut->createPtr();

    (*ssBnd) = -1.0; // clear out field

    Lucee::FieldPtr<double> ssp = ssBnd->createPtr();
    Lucee::FieldPtr<double> sspr = ssBnd->createPtr();

// loop, figuring out where and how the boundary edge is oriented
    for (unsigned dir=0; dir<NDIM; ++dir)
    {
// create sequencer to loop over *each* 1D slice in 'dir' direction
      Lucee::RowMajorSequencer<NDIM> seq(localRgn.deflate(dir));

// lower and upper bounds of 1D slice
      int sliceLower = localRgn.getLower(dir)-1;
      int sliceUpper = localRgn.getUpper(dir);

      while (seq.step())
      {
        int idx[NDIM], idxr[NDIM];
        seq.fillWithIndex(idx);
        seq.fillWithIndex(idxr);
// loop over each slice
        for (int i=sliceLower; i<sliceUpper; ++i)
        {
          idx[dir] = i;
          idxr[dir] = i+1;
// set pointers
          inOut->setPtr(iop, idx);
          inOut->setPtr(iopr, idxr);

          ssBnd->setPtr(ssp, idx);
          ssBnd->setPtr(sspr, idxr);

// flags SSB_UP and SSB_LO indicate direction of ghost cell wrt to a
// skin cell. Note that a cell can "ghost" for left/right AND
// top/bottom AND up/down cells. CAUTION: The degenerate case in which
// a ghost cell is shared by two cells in the same sweep direction is
// not handled. This can occur if the boundary has a very sharp edge,
// like a sharp wedge like object. This degenerate can be handled by
// storing NDIM*2 components in the ssBnd field, but for now I am not
// doing it. (AHH)
// 
//
// (This explanation seems confusing, and should be improved. Ammar
// Hakim, July 8th 2014)

          if ((iop[0] == 1) and (iopr[0] == 0))
// going from inside to outside (edge is on right of cell with index idx)
            ssp[dir] = SSB_UP;
          else if ((iop[0] == 0) and (iopr[0] == 1))
// going from outside to inside (edge is on left of cell with index idxr)
            sspr[dir] = SSB_LO;
        }
      }
    }
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  StairSteppedBcUpdater<NDIM>::update(double t) 
  {
// don't do anything if apply direction is not in range
    if (bcDir>=NDIM)
    {
      Lucee::Except lce("StairSteppedBcUpdater::update: Update direction should be less than ");
      lce << NDIM << ". Provided " << bcDir << " instead";
      throw lce;
    }

    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
// get field to update
    Lucee::Field<NDIM, double>& A = this->getOut<Lucee::Field<NDIM, double> >(0);

// create coordinate system along this direction (note that this
// updater only works on Cartesian meshes. Hence, using an aligned CS
// is perfectly fine)
    Lucee::AlignedRectCoordSys coordSys(bcDir);

    Lucee::FieldPtr<double> Aptr = A.createPtr();
    Lucee::FieldPtr<double> Aptrm = A.createPtr();
    Lucee::ConstFieldPtr<double> ssp = ssBnd->createConstPtr();
    double xc[3];

// create sequencer to loop over *each* 1D slice in 'dir' direction
    Lucee::RowMajorSequencer<NDIM> seq(localRgn.deflate(bcDir));

// lower and upper bounds of 1D slice
    int sliceLower = localRgn.getLower(bcDir);
    int sliceUpper = localRgn.getUpper(bcDir);

    while (seq.step())
    {
      int idx[NDIM], idxm[NDIM];
      seq.fillWithIndex(idx);
      seq.fillWithIndex(idxm);
// loop over each slice
      for (int i=sliceLower; i<sliceUpper; ++i)
      {
        idx[bcDir] = i;
// get centroid coordinate
        grid.setIndex(idx);
        grid.getCentroid(xc);
// set pointers
        A.setPtr(Aptr, idx);
        ssBnd->setPtr(ssp, idx);
          
        if (ssp[bcDir] != -1)
        { // we have found a skin cell
          idxm[bcDir] = i-1;
          if (ssp[bcDir] == SSB_UP)
            idxm[bcDir] = i+1;
// set pointer to right/left cell and apply BC
          A.setPtr(Aptrm, idxm);
          for (std::vector<Lucee::BoundaryCondition*>::const_iterator bcItr = bcList.begin();
               bcItr != bcList.end(); ++bcItr)
            (*bcItr)->applyBc(t, xc, idxm, coordSys, Aptr, Aptrm);
        }
      }
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  StairSteppedBcUpdater<NDIM>::declareTypes() 
  {
// any number of output fields
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void
  StairSteppedBcUpdater<NDIM>::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
// call base class Lua methods
    UpdaterIfc::appendLuaCallableMethods(lfm);

    lfm.appendFunc("setDir", luaSetDir);
  }

  template <unsigned NDIM>
  int
  StairSteppedBcUpdater<NDIM>::luaSetDir(lua_State *L)
  {
    StairSteppedBcUpdater<NDIM> *updater
      = Lucee::PointerHolder<StairSteppedBcUpdater<NDIM> >::getObj(L);
    int d = (unsigned) lua_tonumber(L, 2); // current time to set
    updater->setDir(d);

    return 0;
  }

// instantiations
  template class StairSteppedBcUpdater<1>;
  template class StairSteppedBcUpdater<2>;
  template class StairSteppedBcUpdater<3>;
}
