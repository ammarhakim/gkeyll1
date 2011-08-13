/**
 * @file LcMappedCartGrid.cpp
 *
 * @brief A logically rectangular grid, but non-rectangular in physical space.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcMappedCartGrid.h>
#include <LcStructGridField.h>

namespace Lucee
{
// set ids for grid creators
  template <> const char *MappedCartGrid<1>::id = "MappedCart1D";
  template <> const char *MappedCartGrid<2>::id = "MappedCart2D";
  template <> const char *MappedCartGrid<3>::id = "MappedCart3D";

  template <unsigned NDIM>
  MappedCartGrid<NDIM>::MappedCartGrid()
    : idxr(&Lucee::FixedVector<NDIM, unsigned>(1)[0], &Lucee::FixedVector<NDIM, int>(1)[0])
  {
  }

  template <unsigned NDIM>
  void
  MappedCartGrid<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// get cells in domain
    cells = tbl.getNumVec("cells");
    if (cells.size() != NDIM)
    {
      Lucee::Except lce("MappedCartGrid::readInput: 'cells' should have exactly ");
      lce << NDIM << " elements. Instead has " << cells.size() << std::endl;
      throw lce;
    }

    int ilo[NDIM], iup[NDIM];
    double xlo[NDIM], xup[NDIM];
    for (unsigned i=0; i<NDIM; ++i)
    {
      ilo[i] = 0;
      iup[i] = cells[i];
      xlo[i] = 0.0; // assumes computation region is a unit box
      xup[i] = 1.0;
    }
    Lucee::Region<NDIM, int> localBox(ilo, iup);
    Lucee::Region<NDIM, int> globalBox(ilo, iup);
    Lucee::Region<NDIM, double> physBox(xlo, xup);
// set grid data
    this->setGridData(localBox, globalBox, physBox);

// get vertices array
    Lucee::StructGridField<NDIM, double>& vertices =
      tbl.template getObjectAsDerived<Lucee::StructGridField<NDIM, double> >("vertices");
// get local extended region and indexer
    localExtBox = vertices.getExtRegion();

    idxr = Lucee::RowMajorIndexer<NDIM>(localExtBox);

// allocate vertex data
    unsigned vol = localExtBox.getVolume();
    geometry.setNumVertices(vol);

// copy over data into indexer
    int idx[NDIM];
    Lucee::ConstFieldPtr<double> vPtr = vertices.createConstPtr();
    Lucee::RowMajorSequencer<NDIM> seq(localExtBox);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      vertices.setPtr(vPtr, idx);

      int linIdx = idxr.getIndex(idx);
      for (unsigned k=0; k<NDIM; ++k)
        geometry.vcoords[NDIM*linIdx+k] = vPtr[k];
    }

    Lucee::FixedVector<NDIM, int> zeros(0), minusOnes(-1);
// this actually reduces the size of region by shaving off one layer on upper edges
    Lucee::Region<NDIM, int> geomRegion = localExtBox.extend(&zeros[0], &minusOnes[0]);
    Lucee::RowMajorSequencer<NDIM> geomSeq(geomRegion);

// compute cell geometry
    geometry.setNumCells(vol);
    while (geomSeq.step())
    {
      geomSeq.fillWithIndex(idx);
      int linIdx = idxr.getIndex(idx);
// compute cell geometry based on grid dimension
      if (NDIM == 1)
      {
      }
      else if (NDIM == 2)
      {
      }
      else if (NDIM == 3)
      {
      }
    }

// compute face geometry
    geometry.setNumFaces(vol, true, true);
  }

  template <unsigned NDIM>
  void
  MappedCartGrid<NDIM>::getCentriod(double xc[3]) const
  {
  }

  template <unsigned NDIM>
  void
  MappedCartGrid<NDIM>::getVertex(double xc[3]) const
  {
  }

  template <unsigned NDIM>
  double
  MappedCartGrid<NDIM>::getVolume() const
  {
    return 0.0;
  }

  template <unsigned NDIM>
  double
  MappedCartGrid<NDIM>::getSurfArea(unsigned dir) const
  {
    return 0.0;
  }

  template <unsigned NDIM>
  void
  MappedCartGrid<NDIM>::getSurfCoordSys(unsigned dir, double norm[3],
    double tan1[3], double tan2[3]) const
  {
  }

  template <unsigned NDIM>
  Lucee::IoNodeType
  MappedCartGrid<NDIM>::writeToFile(Lucee::IoBase& io, Lucee::IoNodeType& node,
    const std::string& nm)
  {
// create local and global regions
    Lucee::FixedVector<NDIM, int> zeros(0), ones(1);
    Lucee::Region<NDIM, int> localWriteBox(this->localBox.extend(&zeros[0], &ones[0]));
    Lucee::Region<NDIM, int> globalWriteBox(this->globalBox.extend(&zeros[0], &ones[0]));

// create memory space to write data and copy vertex coordinates (this
// copy is needed (rather than just geometry.vcoords) as
// geometry.vcoords vector has vertex coordinates also for ghost
// cells, which we do not want to write out).
    std::vector<double> buff(NDIM*localWriteBox.getVolume());
    Lucee::RowMajorSequencer<NDIM> seq(localWriteBox);
    unsigned count = 0;
    int idx[NDIM];
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      int linIdx = idxr.getIndex(idx);
      for (unsigned k=0; k<NDIM; ++k)
        buff[count++] = geometry.vcoords[NDIM*linIdx+k];
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
    Lucee::IoNodeType dn =
      io.writeDataSet(node, nm, dataSetSize, dataSetBeg, dataSetLen, &buff[0]);

    io.writeStrAttribute(dn, "vsType", "mesh");
    io.writeStrAttribute(dn, "vsKind", "structured");

    return dn;
  }

  template <unsigned NDIM>
  void
  MappedCartGrid<NDIM>::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
// call base class to register its methods
    Lucee::StructuredGridBase<NDIM>::appendLuaCallableMethods(lfm);
  }

// instantiations
  template class MappedCartGrid<1>;
  template class MappedCartGrid<2>;
  template class MappedCartGrid<3>;
}
