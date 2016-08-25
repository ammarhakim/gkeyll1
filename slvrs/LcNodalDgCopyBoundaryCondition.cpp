/**
 * @file	LcNodalDgCopyBoundaryCondition.cpp
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
#include <LcNodalDgCopyBoundaryCondition.h>

namespace Lucee
{
  static const unsigned LC_LOWER_EDGE = 0;
  static const unsigned LC_UPPER_EDGE = 1;

// set constructor name
  template <> const char *NodalDgCopyBoundaryCondition<1>::id = "NodalDgCopy1D";
  template <> const char *NodalDgCopyBoundaryCondition<2>::id = "NodalDgCopy2D";
  template <> const char *NodalDgCopyBoundaryCondition<3>::id = "NodalDgCopy3D";
  template <> const char *NodalDgCopyBoundaryCondition<4>::id = "NodalDgCopy4D";
  template <> const char *NodalDgCopyBoundaryCondition<5>::id = "NodalDgCopy5D";

  template <unsigned NDIM>
  void
  NodalDgCopyBoundaryCondition<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    BoundaryCondition::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("NodalDgCopyBoundaryCondition::readInput: Must specify element to use using 'basis'");

// read out factors for multiplications
    if (tbl.hasNumVec("fact"))
    {
      fact = tbl.getNumVec("fact");
      if (fact.size() != this->numComponents())
        throw Lucee::Except(
          "NodalDgCopyBoundaryCondition::readInput: If 'fact' table is specified it must have same size as 'components' table");
    }
    else
    {
      fact.resize(this->numComponents());
      for (unsigned i=0; i<fact.size(); ++i) fact[i] = 1.0;
    }
  }

  template <unsigned NDIM>
  void
  NodalDgCopyBoundaryCondition<NDIM>::applyBc(double tm, const double loc[3], const int *idx,
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

    for (unsigned n=0; n<numSurfNodes; ++n)
    {
      int skIdx = skSurfNodes[n]*nc; // start index in skin cell
      int ghIdx = ghSurfNodes[n]*nc; // start index in ghost cell
// copy with scaling factors
      for (unsigned i=0; i<this->numComponents() ; ++i)
        qbc[ghIdx+this->component(i)] = fact[i]*qin[skIdx+this->component(i)];
    }
  }

// instatiations
  template class NodalDgCopyBoundaryCondition<1>;
  template class NodalDgCopyBoundaryCondition<2>;
  template class NodalDgCopyBoundaryCondition<3>;
  template class NodalDgCopyBoundaryCondition<4>;
  template class NodalDgCopyBoundaryCondition<5>;
}
