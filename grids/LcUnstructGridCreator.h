/**
 * @file	LcUnstructGridCreator.h
 *
 * @brief	Unstructured grid creator class.
 *
 * @version	$Id$
 */

#ifndef LC_UNSTRUCT_GRID_CREATOR_H
#define LC_UNSTRUCT_GRID_CREATOR_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// std includes
#include <map>
#include <vector>

// lucee includes
#include <LcUnstructConnectivity.h>
#include <LcUnstructGeometry.h>

namespace Lucee
{
/**
 * Class to create an unstructured grid. This creator can then be
 * passed to UnstructGrid to construct the actual grid. The grid
 * creator needs to be supplied coordinates of each vertex and the
 * cell -> vertex connectivity for each cell. The currently supported
 * cell types are: triangles, quadrilaterals, hexahedrons and
 * tetrahedrons.
 *
 */
  template <typename REAL>
  class UnstructGridCreator
  {
    public:
/**
 * Create a new unstructured grid creator for grid of specified
 * dimension.
 *
 * @param ndim Dimension of grid.
 */
      UnstructGridCreator(unsigned ndim);

/**
 * Set number of vertices in grid.
 *
 * @param nv Number of vertices.
 */
      void setNumVertices(unsigned nv);

/**
 * Set number of cells in grid. A 'cell' is defined as an element of
 * 'ndim' dimension.
 *
 * @param nv Number of vertices.
 */
      void setNumCells(unsigned nc);

/**
 * Append vertex to creator. Vertices are assumed to live in 3D space
 * even for 1D or 2D meshes. Vertices should be added in order that
 * they are numbered.
 *
 * @param xv Coordinates of vertex (x,y,z).
 */
      void appendVertex(double xv[3]);

    private:
/** Dimension of grid */
      unsigned ndim;
/** Vertex coordinates */
      Lucee::UnstructGeometry<3, REAL> vc;
/** Cell->vertex connectivity */
      Lucee::UnstructConnectivity c2v;
/** No copying */
      UnstructGridCreator(const UnstructGridCreator&);
/** No assignment */
      UnstructGridCreator& operator=(const UnstructGridCreator&);
  };
}

#endif // LC_UNSTRUCT_GRID_CREATOR_H
