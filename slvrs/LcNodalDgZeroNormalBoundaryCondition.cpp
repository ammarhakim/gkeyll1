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
      throw Lucee::Except("NodalDisContHyperUpdater::readInput: Must specify element to use using 'basis'");
  }

  template <unsigned NDIM>
  void
  NodalDgZeroNormalBoundaryCondition<NDIM>::applyBc(double tm, const double loc[3], const int *idx,
    const Lucee::RectCoordSys& c, const Lucee::ConstFieldPtr<double>& qin, Lucee::FieldPtr<double>& qbc)
  {
    unsigned numNodes = nodalBasis->getNumNodes();
    unsigned nc = qbc.getNumComponents()/numNodes;

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

    for (unsigned n=0; n<numNodes; ++n)
    {
      // evaluateFunction(*L, tm, nodeCoords, ndIds[n], res);
      // for (unsigned k=0; k<nc; ++k)
      //   qbc[nc*n+k] = res[k];
    }
  }

// instatiations
  template class NodalDgZeroNormalBoundaryCondition<1>;
  template class NodalDgZeroNormalBoundaryCondition<2>;
  template class NodalDgZeroNormalBoundaryCondition<3>;
}
