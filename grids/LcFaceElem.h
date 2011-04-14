/**
 * @file	LcFaceElem.h
 *
 * @brief       Face elements in unstructured grids.
 *
 * @version	$Id$
 */

#ifndef LC_FACE_ELEM_H
#define LC_FACE_ELEM_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// std includes
#include <vector>

namespace Lucee
{
/**
 * A class representing a face.
 */
  template <typename REAL>
  class FaceElem
  {
    public:
/**
 * Fill with face center coordinates.
 *
 * @param xv On output coordinates of face centroid.
 */
      void fillWithCoordinates(REAL xv[3]) const;

/**
 * Get area of face.
 *
 * @return area of face.
 */
      REAL getMeasure() const;

/**
 * Fill with normal (unit vector) to face.
 *
 * @param norm On output normal to face.
 */
      void fillWithNormal(REAL norm[3]) const;

/**
 * Fill with a tangent (unit vector) to face.
 *
 * @param tng On output tangent to face.
 */
      void fillWithTangent1(REAL tng[3]) const;

/**
 * Fill with a tangent (unit) to face. Note that the relation norm =
 * tangent1 X tangent 2 holds.
 *
 * @param tng On output tangent to face.
 */
      void fillWithTangent2(REAL tng[3]) const;

    protected:
/**
 * Create cell given list of cell center coordinates and cell volumes.
 *
 * @param fc Face center coordinates.
 * @param fa Face area.
 */
      FaceElem(const std::vector<REAL>& fc, const std::vector<REAL>& fa);

/**
 * Increment vertex location by one.
 */
      void incr() const
      {
        curr += 1;
      }

/**
 * Are we at end of iteration?
 *
 * @return true
 */
      bool atEnd() const
      {
        if (curr<fa.size())
          return false;
        return true;
      }

    private:
/** Reference to face centroid coordinates */
      const std::vector<REAL>& fc;
/** Reference to face area */
      const std::vector<REAL>& fa;
/** Index of current element */
      mutable unsigned curr;
  };
}

#endif // LC_FACE_ELEM_H
