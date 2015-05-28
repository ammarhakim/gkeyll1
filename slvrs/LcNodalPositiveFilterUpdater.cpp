/**
 * @file	LcNodalPositiveFilterUpdater.cpp
 *
 * @brief	Updater to solve hyperbolic equations with nodal DG scheme.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAlignedRectCoordSys.h>
#include <LcField.h>
#include <LcLinAlgebra.h>
#include <LcMathLib.h>
#include <LcNodalPositiveFilterUpdater.h>
#include <LcStructuredGridBase.h>

// std includes
#include <vector>

namespace Lucee
{
// set id for module system
  template <> const char *NodalPositiveFilterUpdater<1>::id = "NodalPositiveFilter1D";
  template <> const char *NodalPositiveFilterUpdater<2>::id = "NodalPositiveFilter2D";
  template <> const char *NodalPositiveFilterUpdater<3>::id = "NodalPositiveFilter3D";

  template <unsigned NDIM>
  NodalPositiveFilterUpdater<NDIM>::NodalPositiveFilterUpdater()
    : UpdaterIfc()
  {
  }

  template <unsigned NDIM>  
  void 
  NodalPositiveFilterUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("NodalPositiveFilterUpdater::readInput: Must specify element to use using 'basis'");

    if (tbl.hasObject<Lucee::HyperEquation>("equation"))
      equation = &tbl.getObjectAsBase<Lucee::HyperEquation>("equation");
    else
    {
      Lucee::Except lce("NodalPositiveFilterUpdater::readInput: Must specify an equation to solve!");
      throw lce;
    }

    if (tbl.hasString("operation"))
    {
// "flatten" checks cell averages and flattens the cell and its
// neighbors. "filter" checks nodal values and flattens the cell (but
// not its neighbors).
      std::string op = tbl.getString("operation");
      if (op == "flatten")
        opType = OP_FLATTEN;
      else if (op == "filter")
        opType = OP_FILTER;
      else
      {
        Lucee::Except lce("NodalPositiveFilterUpdater::readInput: 'operation' ");
        lce << op << " not recognized!" << std::endl;
        throw lce;
      }
    }
    else
    {
      Lucee::Except lce(
        "NodalPositiveFilterUpdater::readInput: Must specify an operation (one of \"filter\" or \" flatten\") ");
      throw lce;
    }

    nlocal = nodalBasis->getNumNodes();
    meqn = equation->getNumEqns();
  }

  template <unsigned NDIM>
  void 
  NodalPositiveFilterUpdater<NDIM>::initialize()
  {
    Lucee::UpdaterIfc::initialize();
// get weights for quadrature
    weights.resize(nodalBasis->getNumNodes());
    nodalBasis->getWeights(weights);
// get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    double vol = grid.getVolume();
// normalize weights as we are computing averages and not integrals over cell
    for (unsigned i=0; i<nodalBasis->getNumNodes(); ++i)
      weights[i] = weights[i]/vol;
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus 
  NodalPositiveFilterUpdater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

// This may be a bit confusing: qNew is the updated DG solution, which
// is what we check for positivity violation. qOld is the old solution
// which we want to fix (if needed) if positivity is violated.
    const Lucee::Field<NDIM, double>& qNew = this->getInp<Lucee::Field<NDIM, double> >(0);
    Lucee::Field<NDIM, double>& qOld = this->getOut<Lucee::Field<NDIM, double> >(0);

    Lucee::FieldPtr<double> qOldPtr = qOld.createPtr();
    Lucee::ConstFieldPtr<double> qNewPtr = qNew.createConstPtr();

    std::vector<double> qAvg(meqn);
    int idx[NDIM], idx1[NDIM];

    bool isPositive = true; // positive unless otherwise found
    Lucee::Region<NDIM, int> localRgn = qNew.getRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      qNew.setPtr(qNewPtr, idx);

      if (opType == OP_FLATTEN)
      {
// compute average from nodal values
        calcAverage(qNewPtr, qAvg);
// check if positivity is violated
        bool isPos = equation->isInvariantDomain(&qAvg[0]);

        if (isPos == false)
        { 
          isPositive = false;
// positivity is violated, flatten cell ..
          qOld.setPtr(qOldPtr, idx);
          calcAverage(qOldPtr, qAvg);
          resetAllNodes(qOldPtr, qAvg);

// ..  and its neighbors. First, reset "right" cell ...
          for (unsigned d=0; d<NDIM; ++d)
          {
            seq.fillWithIndex(idx1);
            idx1[d] = idx1[d]+1; // right cell
            qOld.setPtr(qOldPtr, idx1);
            calcAverage(qOldPtr, qAvg);
            resetAllNodes(qOldPtr, qAvg);
          }
// ... and then reset "left" cell.
          for (unsigned d=0; d<NDIM; ++d)
          {
            seq.fillWithIndex(idx1);
            idx1[d] = idx1[d]-1;
            qOld.setPtr(qOldPtr, idx1);
            calcAverage(qOldPtr, qAvg);
            resetAllNodes(qOldPtr, qAvg);
          }
        }
      }
      else
      {
        bool isPos = true;
// check if positivity at any node is violated
        for (unsigned k=0; k<nlocal; ++k)
        {
          if (equation->isInvariantDomain(&qNewPtr[k*meqn]) == false)
          { 
            isPos = false; 
            break; 
          }
        }
        if (isPos == false)
        {
          isPositive = false;
// reset this cell
          calcAverage(qOldPtr, qAvg);
          resetAllNodes(qOldPtr, qAvg);     
        }
      }
    }

    return Lucee::UpdaterStatus(isPositive, std::numeric_limits<double>::max());
  }

  template <unsigned NDIM>  
  void
  NodalPositiveFilterUpdater<NDIM>::declareTypes()
  {
// takes one input
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
// returns one output
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void
  NodalPositiveFilterUpdater<NDIM>::calcAverage(const Lucee::ConstFieldPtr<double>& qIn, std::vector<double>& qAvg)
  {
    for (unsigned m=0; m<meqn; ++m)
    {
      qAvg[m] = 0.0;
      for (unsigned k=0; k<nlocal; ++k)
        qAvg[m] += weights[k]*qIn[k*meqn+m];
    }
  }

  template <unsigned NDIM>
  void
  NodalPositiveFilterUpdater<NDIM>::resetAllNodes(Lucee::FieldPtr<double>& q, const std::vector<double>& qAvg)
  {
    for (unsigned k=0; k<nlocal; ++k)
    {
      for (unsigned m=0; m<meqn; ++m)
        q[k*meqn+m] = qAvg[m];
    } 
  }

// instantiations
  template class NodalPositiveFilterUpdater<1>;
  template class NodalPositiveFilterUpdater<2>;
  template class NodalPositiveFilterUpdater<3>;
}
