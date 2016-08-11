/**
 * @file	LcSOLZeroNormalBoundaryCondition.cpp
 *
 * @brief	Class for applying zero gradient boundary conditions on the distribution function
 * Started with LcNodalDgZeroNormalBoundaryCondition and removed coordinate system complications
 *
 * @author eshi
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLZeroNormalBoundaryCondition.h>

namespace Lucee
{
  static const unsigned LC_LOWER_EDGE = 0;
  static const unsigned LC_UPPER_EDGE = 1;

// set constructor name
  template <> const char *SOLZeroNormalBoundaryCondition<5>::id = "SOLZeroNormal5D";

  template <unsigned NDIM>
  void
  SOLZeroNormalBoundaryCondition<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    BoundaryCondition::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("SOLZeroNormalBoundaryCondition::readInput: Must specify element to use using 'basis'");
  }

  template <unsigned NDIM>
  void
  SOLZeroNormalBoundaryCondition<NDIM>::applyBc(double tm, const double loc[3], const int *idx,
    const Lucee::RectCoordSys& c, const Lucee::ConstFieldPtr<double>& qin1,
    const Lucee::ConstFieldPtr<double>& qin, Lucee::FieldPtr<double>& qbc)
  {
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

    // Copy skin cell surface nodes to ghost cell
    for (int nodeIndex = 0; nodeIndex < numSurfNodes; nodeIndex++)
      qbc[ghSurfNodes[nodeIndex]] = qin[skSurfNodes[nodeIndex]];
  }

// instatiations
  template class SOLZeroNormalBoundaryCondition<5>;
}
