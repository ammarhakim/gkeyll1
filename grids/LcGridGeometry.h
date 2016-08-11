/**
 * @file	LcGridGeometry.h
 *
 * @brief	Class holding geometry of unstructured grids.
 */

#ifndef LC_UNSTRUCT_GEOMETRY_H
#define LC_UNSTRUCT_GEOMETRY_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// std includes
#include <vector>

namespace Lucee
{
/**
 * Class to hold unstructured grid geometry information.
 */
  template <unsigned NDIM, typename REAL>
  class GridGeometry
  {
    public:
/**
 * Create empty geometery object. The 'reset' method must be called to
 * allocate any memory to store the geometry.
 */
      GridGeometry();

/**
 * Set number of vertices. This needs to be called before vertex
 * geometry information is set.
 *
 * @param nv Number of vertices.
 */
      void setNumVertices(unsigned nv);

/**
 * Set number of edge. This needs to be called before edge geometry
 * information is set.
 *
 * @param nf Number of edged.
 * @param storeNormal True if normal should be stored.
 * @param storeTangents True if tangents should be stored.
 */
    void setNumEdges(unsigned nf, bool storeNormal, bool storeTangents);

/**
 * Set number of faces. This needs to be called before face geometry
 * information is set.
 *
 * @param nf Number of cells.
 * @param storeNormal True if normal should be stored.
 * @param storeTangents True if tangents should be stored.
 */
      void setNumFaces(unsigned nf, bool storeNormal, bool storeTangents);

/**
 * Set number of cells. This needs to be called before cell geometry
 * information is set.
 *
 * @param nc Number of cells.
 */
      void setNumCells(unsigned nc);

/** Vertex coordinates NDIM*numVertices stored as x,y,z */
      std::vector<REAL> vcoords;

/** Edge center coordinates stored as x,y,z */
      std::vector<REAL> edgeCenter;
/** Edge length */
      std::vector<REAL> edgeLength;
/** Edge normal (only makes sense in 2D) */
      std::vector<REAL> edgeNormal;
/** First edge tangent (only makes sense in 2D) */
      std::vector<REAL> edgeTangent1;
/** Second edge tangent (only makes sense in 2D) */
      std::vector<REAL> edgeTangent2;

/** Face centroid coordinates stored as x,y,z */
      std::vector<REAL> faceCenter;
/** Face area */
      std::vector<REAL> faceArea;
/** Face normal (only makes sense in 3D) */
      std::vector<REAL> faceNormal;
/** First face tangent (only makes sense in 3D) */
      std::vector<REAL> faceTangent1;
/** Second face tangent (only makes sense in 3D) */
      std::vector<REAL> faceTangent2;

/** Cell centroid coordinates stored as x,y,z */
      std::vector<REAL> cellCentroid;
/** Cell volume */
      std::vector<REAL> cellVolume;
  };
}

#endif //  LC_UNSTRUCT_GEOMETRY_H
