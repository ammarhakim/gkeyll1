/**
 * @file	LcUnstructGeometry.h
 *
 * @brief	Class holding geometry of unstructured grids.
 *
 * @version	$Id$
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
  class UnstructGeometry
  {
    public:
/**
 * Create empty geometery object. The 'reset' method must be called to
 * allocate any memory to store the geometry.
 */
      UnstructGeometry();

/**
 * Set number of vertices. This needs to be called before vertex
 * geometry information is set.
 *
 * @param nv Number of vertices.
 */
      void setNumVertices(unsigned nv);

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
/** Face normal */
      std::vector<REAL> faceNormal;
/** First face tangent */
      std::vector<REAL> faceTangent1;
/** Second face tangent */
      std::vector<REAL> faceTangent2;

/** Cell centroid coordinates stored as x,y,z */
      std::vector<REAL> cellCentroid;
/** Cell volume */
      std::vector<REAL> cellVolume;
  };
}

#endif //  LC_UNSTRUCT_GEOMETRY_H
