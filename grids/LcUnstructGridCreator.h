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
// some definitions for cell types
  const short TET_CELL_T=0;
  const short HEX_CELL_T=1;
  const short TRI_CELL_T=2;
  const short QUAD_CELL_T=3;

/**
 * Class to create an unstructured grid. This creator can then be
 * passed to UnstructGrid to construct the actual grid. The grid
 * creator needs to be supplied coordinates of each vertex and the
 * cell -> vertex connectivity for each cell. Here, cell is defined as
 * the element of highest dimension in the mesh. The currently
 * supported cell types are: triangle, quadrilateral, hexahedron
 * and tetrahedron. 
 *
 * Vertices and cells must be numbered starting from 0. To ensure
 * exeception free use, the setNumVertices() and setNumCells() methods
 * must be called first before adding any vertices or cells. Once
 * these methods are called, the vertices can be added in any
 * order. However, the cells must be added one after the other
 * sequentially.
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
 * Get grid dimension.
 *
 * @return grid dimension.
 */
      unsigned getDim() const { return ndim; }

/**
 * Fill with grid vertex information.
 *
 * @param geo Geometry to fill data in.
 */
      void fillWithGeometry(Lucee::UnstructGeometry<3, REAL>& geo) const;

/**
 * Fill with ndim->0 connectivity information.
 *
 * @param conn Connectivity data to fill data in.
 */
      void fillWithConnectivity(Lucee::UnstructConnectivity& conn) const;

/**
 * Get number of triangles in grid.
 *
 * @return number of triangles.
 */
      unsigned getNumTri() const;

/**
 * Get number of quadrilaterals in grid.
 *
 * @return number of quadrilaterals.
 */
      unsigned getNumQuad() const;

/**
 * Get number of tetrahedra in grid.
 *
 * @return number of tetrahedra.
 */
      unsigned getNumTet() const;

/**
 * Get number of hexahedra in grid.
 *
 * @return number of hexahedra.
 */
      unsigned getNumHex() const;

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
 * Set a vertex coordinate. Vertices are assumed to live in 3D space
 * even for 1D or 2D meshes. Vertices should be numbered from
 * 0,1,...,nv-1, where nv is the number of vertices in the grid.
 *
 * @param iv Vertex number.
 * @param xv Coordinates of vertex (x,y,z).
 */
      void setVertex(unsigned iv, REAL xv[3]);

/**
 * Set a vertex coordinate. Vertices are assumed to live in 3D space
 * even for 1D or 2D meshes. Vertices should be numbered from
 * 0,1,...,nv-1, where nv is the number of vertices in the grid.
 *
 * @param iv Vertex number.
 * @param x X-Coordinate of vertex.
 */
      void setVertexXCoord(unsigned iv, REAL x);

/**
 * Set a vertex coordinate. Vertices are assumed to live in 3D space
 * even for 1D or 2D meshes. Vertices should be numbered from
 * 0,1,...,nv-1, where nv is the number of vertices in the grid.
 *
 * @param iv Vertex number.
 * @param y Y-Coordinate of vertex.
 */
      void setVertexYCoord(unsigned iv, REAL y);

/**
 * Set a vertex coordinate. Vertices are assumed to live in 3D space
 * even for 1D or 2D meshes. Vertices should be numbered from
 * 0,1,...,nv-1, where nv is the number of vertices in the grid.
 *
 * @param iv Vertex number.
 * @param z Z-Coordinate of vertex.
 */
      void setVertexZCoord(unsigned iv, REAL z);

/**
 * Add a triangle with vertex numbers (a,b,c). Vertices must be
 * specified in in clockwise order to ensure correct orientation of
 * cells.
 *
 * @param a Vertex index 
 * @param b Vertex index 
 * @param c Vertex index
 */
      void addTri(unsigned a, unsigned b, unsigned c);

/**
 * Add a quad with vertex numbers (a,b,c,d). Vertices must be specified
 * in in clockwise order to ensure correct orientation of cells.
 *
 * @param a Vertex index 
 * @param b Vertex index 
 * @param c Vertex index
 * @param d Vertex index
 */
      void addQuad(unsigned a, unsigned b, unsigned c, unsigned d);

/**
 * Add a tet to the mesh with indices (a,b,c,d). The faces of the tets
 * are [a,b,c],[a,c,d],[b,c,d],[a,b,d]].
 *
 * @param a Vertex index 
 * @param b Vertex index 
 * @param c Vertex index
 * @param d Vertex index
 */
      void addTet(unsigned a, unsigned b, unsigned c, unsigned d);

    private:
/** Dimension of grid */
      unsigned ndim;
/** Current cell number */
      unsigned currCell;
/** number of cells of each type */
      std::map<short, unsigned> cellCount;
/** Vertex coordinates */
      Lucee::UnstructGeometry<3, REAL> vc;
/** Cell->vertex connectivity */
      Lucee::UnstructConnectivity c2v;
/** Cell type */
      std::vector<short> cellType;

/** No copying */
      UnstructGridCreator(const UnstructGridCreator&);
/** No assignment */
      UnstructGridCreator& operator=(const UnstructGridCreator&);
  };
}

#endif // LC_UNSTRUCT_GRID_CREATOR_H
