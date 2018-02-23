/**
 * @file LcNonUniRectCartGrid.cpp
 *
 * @brief A grid allowing for non-uniform spacing, but otherwise rectangular
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGlobals.h>
#include <LcMathLib.h>
#include <LcNonUniRectCartGrid.h>
#include <LcStructGridField.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
// set ids for grid creators
  template <> const char *NonUniRectCartGrid<1>::id = "NonUniformRectCart1D";
  template <> const char *NonUniRectCartGrid<2>::id = "NonUniformRectCart2D";
  template <> const char *NonUniRectCartGrid<3>::id = "NonUniformRectCart3D";
  template <> const char *NonUniRectCartGrid<4>::id = "NonUniformRectCart4D";
  template <> const char *NonUniRectCartGrid<5>::id = "NonUniformRectCart5D";

  template <unsigned NDIM>
  NonUniRectCartGrid<NDIM>::NonUniRectCartGrid()
  {
  }

  template <unsigned NDIM>
  NonUniRectCartGrid<NDIM>::~NonUniRectCartGrid()
  {
    vcoords.clear();
    cellSize.clear();
  }

  template <unsigned NDIM>
  void
  NonUniRectCartGrid<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    StructuredGridBase<NDIM>::readInput(tbl);

    std::vector<int> mapFnRefs;
    std::vector<void *> vcoordsCptrs;
// get list of mapping functions
    if (tbl.hasTable("mappings")) {
      Lucee::LuaTable mapTbl = tbl.getTable("mappings");
      mapFnRefs = mapTbl.getAllFunctionRefs();
    } else if (tbl.hasTable("vertices")) {
      Lucee::LuaTable verticesTbl = tbl.getTable("vertices");
      vcoordsCptrs = verticesTbl.getAllUserdata();
    }

    typename Lucee::Region<NDIM, int> localRgn = this->getLocalRegion();
    typename Lucee::Region<NDIM, double> compSpace = this->getComputationalSpace();

    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;
// compute mappings in each direction
    for (unsigned d=0; d<NDIM; ++d)
    {
// allocate memory to store d-coordinate: we need two ghost cells on
// each side, however, as in a cell we are storing the left vertex
// coordinate, we need one extra cell on the right
      int lower[1], upper[1];
      lower[0] = localRgn.getLower(d)-2;
      upper[0] = localRgn.getUpper(d)+3;
      Lucee::Region<1, int> rgn(lower, upper);
      vcoords.push_back( Array<1, double>(rgn) );

      double low = compSpace.getLower(d);
      double dxc = this->getDx(d); // computational space cell-size

      if (tbl.hasTable("mappings")) {
        for (int i=lower[0]; i<upper[0]; ++i)
        {
          vcoords[d](i) = this->evalLuaFunc(*L, mapFnRefs[d], low+dxc*i);
        }
      } else if (tbl.hasTable("vertices")) {
        double * cptr = (double *) (vcoordsCptrs[d]);
        for (int i=lower[0]; i<upper[0]; ++i)
        {
          vcoords[d](i) = cptr[i+3];
        }
      }
    }

    for (unsigned d=0; d<NDIM; ++d)
    {
      if (this->isPeriodicDir(d))
      {
        // then go fetch the approproate indexes in 
        // GLOBAL regions and reset vertex values.
        int lowerL[1], upperL[1];
        int lowerG[1], upperG[1];
        typename Lucee::Region<NDIM, int> globalRgn = this->getGlobalRegion();
        lowerL[0]  = localRgn.getLower(d);
        upperL[0]  = localRgn.getUpper(d);
        lowerG[0]  = globalRgn.getLower(d);
        upperG[0]  = globalRgn.getUpper(d);
        if (lowerL[0] == lowerG[0]) {
          vcoords[d](lowerL[0]-1) = -(vcoords[d](upperG[0])-vcoords[d](upperG[0]-1));
          vcoords[d](lowerL[0]-2) = -(vcoords[d](upperG[0])-vcoords[d](upperG[0]-2));
        }
        if (upperL[0] == upperG[0]) {
          vcoords[d](upperL[0]+1) = vcoords[d](upperG[0])+(vcoords[d](lowerG[0]+1)-vcoords[d](lowerG[0]));
          vcoords[d](upperL[0]+2) = vcoords[d](upperG[0])+(vcoords[d](lowerG[0]+2)-vcoords[d](lowerG[0]));
        }
      }
    }

// compute cell sizes
    for (unsigned d=0; d<NDIM; ++d)
    {
      int lower[1], upper[1];
      lower[0] = localRgn.getLower(d)-2;
      upper[0] = localRgn.getUpper(d)+2;
      Lucee::Region<1, int> rgn(lower, upper);
      cellSize.push_back( Array<1, double>(rgn) );
      for (int i=lower[0]; i<upper[0]; ++i)
        cellSize[d](i) = vcoords[d](i+1)-vcoords[d](i);
    }
  }

  template <unsigned NDIM>
  void
  NonUniRectCartGrid<NDIM>::getCentroid(double xc[]) const
  {
    for (unsigned d=0; d<NDIM; ++d)
      xc[d] = 0.5*(vcoords[d](this->currIdx[d]+1) + vcoords[d](this->currIdx[d]));
    for (unsigned d=NDIM; d<3; ++d)
      xc[d] = 0.0;
  }

  template <unsigned NDIM>
  void
  NonUniRectCartGrid<NDIM>::getVertex(double xc[]) const
  {
    for (unsigned d=0; d<NDIM; ++d)
      xc[d] = vcoords[d](this->currIdx[d]);
    for (unsigned d=NDIM; d<3; ++d)
      xc[d] = 0.0;
  }

  template <unsigned NDIM>
  double
  NonUniRectCartGrid<NDIM>::getVolume() const
  {
    double vol = 1.0;
    for (unsigned d=0; d<NDIM; ++d)
      vol *= cellSize[d](this->currIdx[d]);
    return vol;
  }

  template <unsigned NDIM>
  double
  NonUniRectCartGrid<NDIM>::getSurfArea(unsigned dir) const
  {
    double sa = 1.0;
    for (unsigned d=0; d<NDIM; ++d)
    {
      if (dir!=d)
        sa *= cellSize[d](this->currIdx[d]);
    }
    return sa;
  }

  template <unsigned NDIM>
  void
  NonUniRectCartGrid<NDIM>::getSurfCoordSys(unsigned dir, double norm[],
    double tan1[], double tan2[]) const
  {
    if (dir==0)
    {
      norm[0] = 1.0;
      norm[1] = 0.0;
      norm[2] = 0.0;

      tan1[0] = 0.0;
      tan1[1] = 1.0;
      tan1[2] = 0.0;

      tan2[0] = 0.0;
      tan2[1] = 0.0;
      tan2[2] = 1.0;
    }
    else if (dir==1)
    {
      norm[0] = 0.0;
      norm[1] = 1.0;
      norm[2] = 0.0;

      tan1[0] = -1.0;
      tan1[1] = 0.0;
      tan1[2] = 0.0;

      tan2[0] = 0.0;
      tan2[1] = 0.0;
      tan2[2] = 1.0;
    }
    else if (dir==2)
    {
      norm[0] = 0.0;
      norm[1] = 0.0;
      norm[2] = 1.0;

      tan1[0] = 1.0;
      tan1[1] = 0.0;
      tan1[2] = 0.0;

      tan2[0] = 0.0;
      tan2[1] = 1.0;
      tan2[2] = 0.0;
    }
  }

  template <unsigned NDIM>
  TxIoNodeType
  NonUniRectCartGrid<NDIM>::writeToFile(TxIoBase& io, TxIoNodeType& node,
    const std::string& nm)
  {
// NOTE: Even though VizSchema allows a non-uniform grid, I am writing
// it out as if it is a general structured (mapped) grid. Makes
// post-processing tools easier to maintain at the expense of a
// somewhat larger amount of data writen. (Ammar Hakim, April 21st
// 2014).

// create local and global regions
    Lucee::FixedVector<NDIM, int> zeros(0), ones(1);
    Lucee::Region<NDIM, int> localWriteBox(this->localRgn.extend(&zeros[0], &ones[0]));
    Lucee::Region<NDIM, int> globalWriteBox(this->globalRgn.extend(&zeros[0], &ones[0]));

// create memory space to write data and copy vertex coordinates
    std::vector<double> buff(NDIM*localWriteBox.getVolume());
    Lucee::RowMajorSequencer<NDIM> seq(localWriteBox);
    unsigned count = 0;
    int idx[NDIM];
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      for (unsigned k=0; k<NDIM; ++k)
        buff[count++] = vcoords[k](idx[k]);
    }

    std::vector<size_t> dataSetSize, dataSetBeg, dataSetLen;
// construct sizes and shapes to write stuff out
    for (unsigned i=0; i<NDIM; ++i)
    {
      dataSetSize.push_back( globalWriteBox.getShape(i) );
      dataSetBeg.push_back( localWriteBox.getLower(i) - globalWriteBox.getLower(i) );
      dataSetLen.push_back( localWriteBox.getShape(i) );
    }
    dataSetSize.push_back(NDIM);
    dataSetBeg.push_back(0);
    dataSetLen.push_back(NDIM);

// write it out
    TxIoNodeType dn =
      io.writeDataSet(node, nm, dataSetSize, dataSetBeg, dataSetLen, &buff[0]);

    io.writeAttribute(dn, "vsType", "mesh");
    io.writeAttribute(dn, "vsKind", "structured");

    // Code copied from LcRectCartGrid to make output fields
    // compatible with RectCart plotting routines
    std::vector<double> lower(NDIM), upper(NDIM);
    std::vector<unsigned> numPhysCells(NDIM), start(NDIM);
    for (unsigned i=0; i<NDIM; ++i)
    {
      lower[i] = this->compSpace.getLower(i);
      upper[i] = this->compSpace.getUpper(i);
      start[i] = this->globalRgn.getLower(i);
      numPhysCells[i] = this->globalRgn.getShape(i);
    }

    io.template 
      writeAttribute<unsigned>(dn, "vsNumCells", numPhysCells);
    io.template
      writeAttribute<double>(dn, "vsLowerBounds", lower);
    io.template
      writeAttribute<double>(dn, "vsUpperBounds", upper);

    return dn;
  }

  template <unsigned NDIM>
  void
  NonUniRectCartGrid<NDIM>::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
// call base class to register its methods
    Lucee::StructuredGridBase<NDIM>::appendLuaCallableMethods(lfm);
  }

  template <unsigned NDIM>
  double
  NonUniRectCartGrid<NDIM>::evalLuaFunc(Lucee::LuaState& L, int fnRef, double inp)
  {
// push function object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, fnRef);
// push variables on stack
    lua_pushnumber(L, inp);
    if (lua_pcall(L, 1, 1, 0) != 0)
    {
      Lucee::Except lce("NonUniRectCartGrid: ");
      lce << "Problem evaluating function supplied to 'mappings' function list";
      throw lce;
    }
// fetch results
    if (!lua_isnumber(L, -1))
      throw Lucee::Except("NonUniRectCartGrid: Return value not a number");
    double res = lua_tonumber(L, -1);
    lua_pop(L, 1);
    return res;
  }

// instantiations
  template class NonUniRectCartGrid<1>;
  template class NonUniRectCartGrid<2>;
  template class NonUniRectCartGrid<3>;
  template class NonUniRectCartGrid<4>;
  template class NonUniRectCartGrid<5>;
}
