/**
 * @file	LcNodalDgZeroNormalBoundaryCondition.cpp
 *
 * @brief	Class for applying copy boundary conditions.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGlobals.h>
#include <LcLuaState.h>
#include <LcNodalDgZeroNormalBoundaryCondition.h>

namespace Lucee
{
  static const unsigned LC_LOWER_EDGE = 0;
  static const unsigned LC_UPPER_EDGE = 1;

// set constructor name
  template <> const char *NodalDgZeroNormalBoundaryCondition<1>::id = "NodalDgZeroNormal1D";
  template <> const char *NodalDgZeroNormalBoundaryCondition<2>::id = "NodalDgZeroNormal2D";
  template <> const char *NodalDgZeroNormalBoundaryCondition<3>::id = "NodalDgZeroNormal3D";

  template <unsigned NDIM>
  void
  NodalDgZeroNormalBoundaryCondition<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    BoundaryCondition::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("NodalDgZeroNormalBoundaryCondition::readInput: Must specify element to use using 'basis'");
  }

  template <unsigned NDIM>
  void
  NodalDgZeroNormalBoundaryCondition<NDIM>::applyBc(double tm, const double loc[3], const int *idx,
    const Lucee::RectCoordSys& c, const Lucee::ConstFieldPtr<double>& qin1,
    const Lucee::ConstFieldPtr<double>& qin, Lucee::FieldPtr<double>& qbc)
  {
    unsigned numNodes = nodalBasis->getNumNodes();
    unsigned nc = qin.getNumComponents()/numNodes;

    unsigned numSurfNodes = nodalBasis->getNumSurfLowerNodes(this->getDir());

    std::vector<int> skSurfNodes(numSurfNodes), ghSurfNodes(numSurfNodes);
// get surface nodes for skin/ghost cells
    nodalBasis->setIndex(idx);
    if (this->getEdge() == LC_LOWER_EDGE)
    {
// lower nodes in skin copied to upper nodes in ghost
      nodalBasis->getSurfLowerNodeNums(this->getDir(), skSurfNodes);
      nodalBasis->getSurfUpperNodeNums(this->getDir(), ghSurfNodes);
    }
    else
    {
// upper nodes in skin copied to lowr nodes in ghost
      nodalBasis->getSurfUpperNodeNums(this->getDir(), skSurfNodes);
      nodalBasis->getSurfLowerNodeNums(this->getDir(), ghSurfNodes);
    }

    double vel[3], normVel[3];
    for (unsigned n=0; n<numSurfNodes; ++n)
    {
      int skIdx = skSurfNodes[n]*nc; // start index in skin cell
// rotate vector to local coordinate system
      for (unsigned i=0; i<3; ++i)
        vel[i] = qin[skIdx+this->component(i)];
      c.rotateVecToLocal(vel, normVel);
      
// flip sign of normal component
      normVel[0] = -normVel[0];

      int ghIdx = ghSurfNodes[n]*nc; // start index in ghost cell
// rotate back to global frame
      c.rotateVecToGlobal(normVel, vel);
        for (unsigned i=0; i<3; ++i)
          qbc[ghIdx+this->component(i)] = vel[i];
    }
  }

// instatiations
  template class NodalDgZeroNormalBoundaryCondition<1>;
  template class NodalDgZeroNormalBoundaryCondition<2>;
  template class NodalDgZeroNormalBoundaryCondition<3>;
}
