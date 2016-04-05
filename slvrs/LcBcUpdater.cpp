/**
 * @file	LcBcUpdater.cpp
 *
 * @brief	Apply boundary conditions.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAlignedRectCoordSys.h>
#include <LcBcUpdater.h>
#include <LcField.h>
#include <LcStructuredGridBase.h>

namespace Lucee
{
  static const unsigned LC_LOWER_EDGE = 0;
  static const unsigned LC_UPPER_EDGE = 1;

  template <> const char *BcUpdater<1>::id = "Bc1D";
  template <> const char *BcUpdater<2>::id = "Bc2D";
  template <> const char *BcUpdater<3>::id = "Bc3D";
  template <> const char *BcUpdater<4>::id = "Bc4D";
  template <> const char *BcUpdater<5>::id = "Bc5D";

  template <unsigned NDIM>
  BcUpdater<NDIM>::BcUpdater() 
    : Lucee::UpdaterIfc() 
  {
  }

  template <unsigned NDIM>
  BcUpdater<NDIM>::~BcUpdater() 
  {
  }

  template <unsigned NDIM>
  void
  BcUpdater<NDIM>::readInput(Lucee::LuaTable& tbl) 
  {
    UpdaterIfc::readInput(tbl);
// read direction to apply
    dir = tbl.getNumber("dir");
    std::string edgeStr = tbl.getString("edge");
    if (edgeStr == "lower")
      edge = LC_LOWER_EDGE;
    else
      edge = LC_UPPER_EDGE;

// get list of boundary conditions to apply
    if ( !tbl.hasTable("boundaryConditions") )
      Lucee::Except lce("BcUpdater::readInput: Must specify boundary conditions to apply!");

    Lucee::LuaTable bcTbl = tbl.getTable("boundaryConditions");
    bcList = bcTbl.template getAllObjects<Lucee::BoundaryCondition>();
  }

  template <unsigned NDIM>
  void
  BcUpdater<NDIM>::initialize() 
  {
// call base class method
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  BcUpdater<NDIM>::update(double t) 
  {
    int lo[NDIM], up[NDIM];

    const Lucee::StructuredGridBase<NDIM>& grid
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
// create coordinate system along this direction
    Lucee::AlignedRectCoordSys coordSys(dir);

    int idx[NDIM], idxG[NDIM];
    double xc[3];
// loop over each array and apply boundary conditions
    for (unsigned n=0; n<this->getNumOutVars(); ++n) {
// get array
      Lucee::Field<NDIM, double>& A = this->getOut<Lucee::Field<NDIM, double> >(n);
      Lucee::ConstFieldPtr<double> iPtr = A.createConstPtr(); // skin
      Lucee::ConstFieldPtr<double> iPtr1 = A.createConstPtr(); // "left" of skin
      Lucee::FieldPtr<double> gPtr = A.createPtr();

// create a region to represent ghost layer
      for (unsigned i=0; i<NDIM; ++i)
      { // whole region, including extended region
        lo[i] = A.getGlobalLowerExt(i);
        up[i] = A.getGlobalUpperExt(i);
      }
// adjust region so it only indexes ghost cells
      if (edge == LC_LOWER_EDGE)
      { // lower side
        up[dir] = A.getGlobalLower(dir);
      }
      else
      { // upper side
        lo[dir] = A.getGlobalUpper(dir);
      }
// region must be local to processor
      Lucee::Region<NDIM, int> gstRgn = A.getExtRegion().intersect(
        Lucee::Region<NDIM, int>(lo, up));

      Lucee::RowMajorSequencer<NDIM> seq(gstRgn);
      while (seq.step())
      {
        seq.fillWithIndex(idx);
        seq.fillWithIndex(idxG);
// get centroid coordinate
        grid.setIndex(idx);
        grid.getCentroid(xc);

        A.setPtr(gPtr, idx);
// set pointer to skin cell
        if (edge == LC_LOWER_EDGE)
          idx[dir] = A.getLower(dir);
        else
          idx[dir] = A.getUpper(dir)-1;
        A.setPtr(iPtr, idx);
// set pointer to "left" of skin cell
        if (edge == LC_LOWER_EDGE)
          idx[dir] = A.getLower(dir)+1;
        else
          idx[dir] = A.getUpper(dir)-2;
        A.setPtr(iPtr1, idx);
// apply boundary conditions
        for (std::vector<Lucee::BoundaryCondition*>::const_iterator bcItr = bcList.begin();
             bcItr != bcList.end(); ++bcItr)
        {
          (*bcItr)->setDir(dir);
          (*bcItr)->setEdge(edge);
          (*bcItr)->applyBc(t, xc, idxG, coordSys, iPtr1, iPtr, gPtr);
        }
      }
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  BcUpdater<NDIM>::declareTypes() 
  {
// any number of output fields
    this->setLastInpVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class BcUpdater<1>;
  template class BcUpdater<2>;
  template class BcUpdater<3>;
  template class BcUpdater<4>;
  template class BcUpdater<5>;
}
