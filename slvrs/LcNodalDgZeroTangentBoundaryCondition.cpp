/**
 * @file	LcNodalDgZeroTangentBoundaryCondition.cpp
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
#include <LcNodalDgZeroTangentBoundaryCondition.h>

namespace Lucee
{
  static const unsigned LC_LOWER_EDGE = 0;
  static const unsigned LC_UPPER_EDGE = 1;

// set constructor name
  template <> const char *NodalDgZeroTangentBoundaryCondition<1>::id = "NodalDgZeroTangent1D";
  template <> const char *NodalDgZeroTangentBoundaryCondition<2>::id = "NodalDgZeroTangent2D";
  template <> const char *NodalDgZeroTangentBoundaryCondition<3>::id = "NodalDgZeroTangent3D";

  template <unsigned NDIM>
  void
  NodalDgZeroTangentBoundaryCondition<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    BoundaryCondition::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("NodalDgZeroTangentBoundaryCondition::readInput: Must specify element to use using 'basis'");
  }

  template <unsigned NDIM>
  void
  NodalDgZeroTangentBoundaryCondition<NDIM>::applyBc(double tm, const double loc[3], const int *idx,
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
      
// flip sign of tangent components
    normVel[1] = -normVel[1];
    normVel[2] = -normVel[2];

      int ghIdx = ghSurfNodes[n]*nc; // start index in ghost cell
// rotate back to global frame
      c.rotateVecToGlobal(normVel, vel);
        for (unsigned i=0; i<3; ++i)
          qbc[ghIdx+this->component(i)] = vel[i];
    }
  }

// instatiations
  template class NodalDgZeroTangentBoundaryCondition<1>;
  template class NodalDgZeroTangentBoundaryCondition<2>;
  template class NodalDgZeroTangentBoundaryCondition<3>;
}
