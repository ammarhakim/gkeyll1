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
#include <LcMathLib.h>
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
// call base class method
    StructuredGridBase<NDIM>::readInput(tbl);

// get vertices array
    Lucee::StructGridField<NDIM, double>& vertices =
      tbl.template getObjectAsDerived<Lucee::StructGridField<NDIM, double> >("vertices");
// get local extended region and indexer
    localExtBox = vertices.getExtRegion();
    Lucee::Region<NDIM, int> localVBox = vertices.getRegion();
// ensure there are correct number of vertices (i.e. vertiex
// coordinates of ghost cells have been specified)
    for (unsigned i=0; i<NDIM; ++i)
    {
      if ((localExtBox.getUpper(i)-localVBox.getUpper(i) != 3) ||
        (localVBox.getLower(i)-localExtBox.getLower(i) != 2))
      {
        throw Lucee::Except("Lucee::MappedCartGrid: 'vertices' field should have {2, 3} ghost cells");
      }
    }

    idxr = Lucee::RowMajorIndexer<NDIM>(localExtBox);

// allocate vertex data
    unsigned vol = localExtBox.getVolume();
    geometry.setNumVertices(vol);

// copy over data into geometry structure
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
// this reduces size of region by shaving off one layer on upper edges
    Lucee::Region<NDIM, int> geomRegion = localExtBox.extend(&zeros[0], &minusOnes[0]);
    Lucee::RowMajorSequencer<NDIM> geomSeq(geomRegion);

// compute grid geometry information
    geometry.setNumCells(vol);
    geometry.setNumFaces(vol, true, true); // always store tangents and normals

    int idx1[NDIM], idx2[NDIM], idx3[NDIM];
    while (geomSeq.step())
    {
      geomSeq.fillWithIndex(idx);
      int linIdx = idxr.getIndex(idx);
// compute cell geometry based on grid dimension
      if (NDIM == 1)
      {
        Lucee::Vec3<double> a(geometry.vcoords[1*linIdx], 0.0, 0.0);
        Lucee::Vec3<double> b(geometry.vcoords[1*linIdx], 0.0, 0.0);
        calc1dGeom(a, b);
      }
      else if (NDIM == 2)
      {
      }
      else if (NDIM == 3)
      {
      }
    }

  }

  template <unsigned NDIM>
  void
  MappedCartGrid<NDIM>::getCentroid(double xc[3]) const
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
  TxIoNodeType
  MappedCartGrid<NDIM>::writeToFile(TxIoBase& io, TxIoNodeType& node,
    const std::string& nm)
  {
// create local and global regions
    Lucee::FixedVector<NDIM, int> zeros(0), ones(1);
    Lucee::Region<NDIM, int> localWriteBox(this->localRgn.extend(&zeros[0], &ones[0]));
    Lucee::Region<NDIM, int> globalWriteBox(this->globalRgn.extend(&zeros[0], &ones[0]));

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
    TxIoNodeType dn =
      io.writeDataSet(node, nm, dataSetSize, dataSetBeg, dataSetLen, &buff[0]);

    io.writeAttribute(dn, "vsType", "mesh");
    io.writeAttribute(dn, "vsKind", "structured");

    return dn;
  }

  template <unsigned NDIM>
  void
  MappedCartGrid<NDIM>::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
// call base class to register its methods
    Lucee::StructuredGridBase<NDIM>::appendLuaCallableMethods(lfm);
  }

  template <unsigned NDIM>
  void
  MappedCartGrid<NDIM>::calc1dGeom(const Lucee::Vec3<double>& a, const Lucee::Vec3<double>& b) const
  {
  }

  template <unsigned NDIM>
  void
  MappedCartGrid<NDIM>::calc2dGeom(const Lucee::Vec3<double>& a, const Lucee::Vec3<double>& b,
    const Lucee::Vec3<double>& c, const Lucee::Vec3<double>& d) const
  {
  }

// instantiations
  template class MappedCartGrid<1>;
  template class MappedCartGrid<2>;
  template class MappedCartGrid<3>;
}
