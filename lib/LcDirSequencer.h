/**
 * @file	LcDirSequencer.h
 *
 * @brief	Sequence over region using row-major order.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_DIR_SEQUENCER_H
#define LC_DIR_SEQUENCER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcExcept.h>
#include <LcRegion.h>

// std includes
#include <vector>

namespace Lucee
{
/**
 * Class to sequence over region, in sweeps in specified direction.
 */
  template <unsigned NDIM>
  class DirSequencer
  {
    public:
/**
 * Create a new directional sequencer.
 *
 * @param rgn Region to sequence.
 * @param dir Direction to sequence in.
 * @param nl Number of cells on the left to include.
 * @param nr Number of cells on the right to include.
 */
      DirSequencer(const Lucee::Region<NDIM, int>& rgn, unsigned dir, unsigned nl, unsigned nr);

/**
 * Moves sequencer by one element. Returns false when there are no
 * more elements.
 */
      bool step();

/**
 * Reset the sequence to run again.
 */
      void reset();

/**
 * Get index at the specified location in stencil.
 *
 * @param loc Location in stencil 
 *    (0 for current, < 0 for left an > 0 for right).
 * @param idx (out) Index at stencil location.
 */
      void fillWithIndex(int loc, int idx[NDIM]) const;

/**
 * Was a stencil location already visited?
 *
 * @param loc Location in stencil 
 *    (0 for current, < 0 for left an > 0 for right).
 * @return true if location visited, false otherwise.
 */
      bool isLocVisited(int loc) const 
      {
        return visited[nl+loc];
      }

    private:
/** Region to sequence */
      Lucee::Region<NDIM, int> rgn;
/** Primary direction */
      unsigned dir;
/** Number of stencil points to the left */
      unsigned nl;
/** Number of stencil points to the right */
      unsigned nr;
/** Flag to indicate if points have been visited */
      std::vector<bool> visited;
/** Current index */
      int index[NDIM];
/** Is this the first time step() is called */
      bool isFirst;
/** Is this box empty? */
      bool isEmpty;
/** Order in which indices must be incremeneted */
      unsigned incOrder[NDIM];

/**
 * Set visited array to false.
 */
      void setVisitedToFalse();
  };

  template <unsigned NDIM>
  DirSequencer<NDIM>::DirSequencer(const Lucee::Region<NDIM, int>& rgn, unsigned dir, unsigned nl, unsigned nr)
    : rgn(rgn), dir(dir), nl(nl), nr(nr), visited(nr+nl+1), isFirst(true), isEmpty(false)
  {
    if (dir>=NDIM)
      throw Lucee::Except("DirSequencer::DirSequencer: direction must be less than dimentions");

// by default no points have been visited
    setVisitedToFalse();
// set initial set of indices
    for (unsigned i=0; i<NDIM; ++i)
      index[i] = rgn.getLower(i);
    if (rgn.getVolume() == 0)
      isEmpty = false; // nothing to do for empty rgn

// construct an array for order in which to increment indices. This is
// essential the cyclic permutation of (0,1,...,NDIM-1) whose first
// element is dir.
    for (unsigned i=0; i<NDIM-dir; ++i)
      incOrder[i] = dir+i;
    for (unsigned i=0; i<dir; ++i)
      incOrder[NDIM-dir+i] = i;
  }

  template <unsigned NDIM>
  void
  DirSequencer<NDIM>::fillWithIndex(int loc, int idx[NDIM]) const {
    for (unsigned i=0; i<NDIM; ++i)
      idx[i] = index[i];
    idx[dir] = index[dir]+loc;
  }

  template <unsigned NDIM>
  bool
  DirSequencer<NDIM>::step() {
    if (isEmpty) return false;

    if (isFirst) {
// first time around: indices are set in ctor
      isFirst = false;
      return true;
    }

    for (unsigned i=0; i<(nr+nl); ++i)
      visited[i] = true;

// simply increment indices till we run out of indices to increment
// and then bump index along next direction.
    unsigned incDir;
    for (unsigned i=0; i<NDIM; ++i) {
      incDir = incOrder[i]; // current direction we are working on
      index[incDir] += 1;
      if (index[incDir] > rgn.getUpper(incDir)-1) {
        setVisitedToFalse();
        index[incDir] = rgn.getLower(incDir);
      }
      else {
        return true;
      }
    }
    return false;
  }

  template <unsigned NDIM>
  void
  DirSequencer<NDIM>::reset() {
    isFirst = true;
    for (unsigned i=0; i<NDIM; ++i)
      index[i] = rgn.getLower(i);
    setVisitedToFalse();
  }

  template <unsigned NDIM>
  void
  DirSequencer<NDIM>::setVisitedToFalse() 
  {
    for (unsigned i=0; i<(nr+nl+1); ++i)
      visited[i] = false;
  }
}

#endif // LC_DIR_SEQUENCER_H
