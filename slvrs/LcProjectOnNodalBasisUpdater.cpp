/**
 * @file	LcProjectOnNodalBasisUpdater.cpp
 *
 * @brief	Project a function of a basis functions.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcProjectOnNodalBasisUpdater.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <vector>

// eigen includes
#include <Eigen/Dense>
#include <Eigen/LU>

namespace Lucee
{

  template <> const char *ProjectOnNodalBasisUpdater<1>::id = "ProjectOnNodalBasis1D";
  template <> const char *ProjectOnNodalBasisUpdater<2>::id = "ProjectOnNodalBasis2D";
  template <> const char *ProjectOnNodalBasisUpdater<3>::id = "ProjectOnNodalBasis3D";
  template <> const char *ProjectOnNodalBasisUpdater<4>::id = "ProjectOnNodalBasis4D";
  template <> const char *ProjectOnNodalBasisUpdater<5>::id = "ProjectOnNodalBasis5D";

  template <unsigned NDIM>
  ProjectOnNodalBasisUpdater<NDIM>::ProjectOnNodalBasisUpdater()
    : fnRef(-1), sharedNodes(false)
  {
  }

  template <unsigned NDIM>
  void
  ProjectOnNodalBasisUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    UpdaterIfc::readInput(tbl);
// get function to evaluate
    fnRef = tbl.getFunctionRef("evaluate");

    sharedNodes = false;
// check if there are shared nodes
    if (tbl.hasBool("shareCommonNodes"))
      sharedNodes = tbl.getBool("shareCommonNodes");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("ProjectOnNodalBasisUpdater::readInput: Must specify element to use using 'basis'");
  }

  template <unsigned NDIM>
  void
  ProjectOnNodalBasisUpdater<NDIM>::initialize()
  {
// call base class method
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  ProjectOnNodalBasisUpdater<NDIM>::update(double t)
  {
// get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
// get output array
    Lucee::Field<NDIM, double>& q = this->getOut<Lucee::Field<NDIM, double> >(0);
    q = 0.0; // set all entries to 0.0

// get list of nodes exclusively owned by element
    std::vector<int> ndIds;

    if (sharedNodes)
      nodalBasis->getExclusiveNodeIndices(ndIds);
    else
    { // create "unit" mapping
      ndIds.resize(nodalBasis->getNumNodes());
      for (unsigned i = 0; i<ndIds.size(); ++i)
        ndIds[i] = i;
    }

// number of nodes
    unsigned numNodes = ndIds.size();

// determine number of components in field
    unsigned nc = q.getNumComponents()/numNodes;

    std::vector<double> res(nc); // to store function evaluation result

// indices into grid
    int idx[NDIM];

// get hold of Lua state object
    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;

// pointer to data in field
    Lucee::FieldPtr<double> ptr = q.createPtr();

// loop over each cell in extended region
    Lucee::Region<NDIM, int> localExtRgn = q.getExtRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localExtRgn);

    // Needed variables in loop:
    int numLocalNodes = nodalBasis->getNumNodes();
    Lucee::Matrix<double> localMass(numLocalNodes,numLocalNodes);
    Eigen::MatrixXd localMassEigen(numLocalNodes,numLocalNodes);
    Eigen::MatrixXd localMassInv(numLocalNodes,numLocalNodes);
    unsigned numVolQuadNodes = nodalBasis->getNumGaussNodes();
    Lucee::Matrix<double> interpMatrix(numVolQuadNodes,numLocalNodes);
    Lucee::Matrix<double> gaussOrdinates(numVolQuadNodes,NC);
    std::vector<double> gaussWeights(numVolQuadNodes);
    Eigen::MatrixXd basisIntegrals(numLocalNodes,nc);
    Eigen::MatrixXd projectF(numLocalNodes,nc);
    // Scaling factors
    double coordScales[NDIM];
    for (unsigned d = 0; d < NDIM; d++)
      coordScales[d] = 0.5*grid.getDx(d);
    double xc[NC];

    while (seq.step())
    {
      seq.fillWithIndex(idx);
      q.setPtr(ptr, idx);

      nodalBasis->setIndex(idx);

      // Get the mass matrix and invert it
      nodalBasis->getMassMatrix(localMass);
      // Copy over to Eigen matrix
      for (int rowIndex = 0; rowIndex < localMassEigen.rows(); rowIndex++)
        for (int colIndex = 0; colIndex < localMassEigen.cols(); colIndex++)
          localMassEigen(rowIndex,colIndex) = localMass(rowIndex,colIndex);

      localMassInv = localMassEigen.inverse();

      // Get interpolation matrix, gaussian quadrature points, and weights
      nodalBasis->getGaussQuadData(interpMatrix,gaussOrdinates,gaussWeights);
      // Translate and scale gaussOrdinates to physical space
      grid.setIndex(idx);
      grid.getCentroid(xc);
      for (int gaussPoint = 0; gaussPoint < numVolQuadNodes; gaussPoint++)
        for (int dimIndex = 0; dimIndex < NC; dimIndex ++)
          gaussOrdinates(gaussPoint,dimIndex) = xc[dimIndex] + gaussOrdinates(gaussPoint,dimIndex)*coordScales[dimIndex];

      basisIntegrals.setZero(basisIntegrals.rows(),basisIntegrals.cols());

      // Evaluate integral of each basis function * function
      // Loop over all basis functions
      for (unsigned gaussPoint = 0; gaussPoint < numVolQuadNodes; gaussPoint++)
      {
        // Get value of function at each gauss quadrature ordinate
        evaluateFunction(*L, t, gaussOrdinates, gaussPoint, res);
        
        for (int basisIndex = 0; basisIndex < basisIntegrals.rows(); basisIndex++)
          for (int component = 0; component < nc; component++)
            basisIntegrals(basisIndex,component) += gaussWeights[gaussPoint]*interpMatrix(gaussPoint,basisIndex)*
              res[component];
      }

      // Columns are components of f (can be more than 1)
      projectF = localMassInv*basisIntegrals;

      for (unsigned nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
      {
        // Copy result into field
        for (unsigned component = 0; component < nc; component++)
          ptr[nc*component+nodeIndex] = projectF(nodeIndex,component);
      }
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  ProjectOnNodalBasisUpdater<NDIM>::declareTypes()
  {
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void
  ProjectOnNodalBasisUpdater<NDIM>::evaluateFunction(Lucee::LuaState& L, double tm,
    const Lucee::Matrix<double> nc, unsigned nn, std::vector<double>& res)
  {
// push function object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, fnRef);
// push variables on stack
    for (unsigned i=0; i<NC; ++i)
      lua_pushnumber(L, nc(nn,i));
    lua_pushnumber(L, tm);
// call function
    if (lua_pcall(L, NC+1, res.size(), 0) != 0)
    {
      Lucee::Except lce("ProjectOnNodalBasisUpdater::evaluateFunction: ");
      lce << "Problem evaluating function supplied as 'evaluate' "
          << std::endl;
      std::string err(lua_tostring(L, -1));
      lua_pop(L, 1);
      lce << "[" << err << "]";
      throw lce;
    }
// fetch results
    for (int i=-res.size(); i<0; ++i)
    {
      if (!lua_isnumber(L, i))
        throw Lucee::Except("ProjectOnNodalBasisUpdater::evaluateFunction: Return value not a number");
      res[res.size()+i] = lua_tonumber(L, i);
    }
    lua_pop(L, 1);
  }

// instantiations
  template class ProjectOnNodalBasisUpdater<1>;
  template class ProjectOnNodalBasisUpdater<2>;
  template class ProjectOnNodalBasisUpdater<3>;
  template class ProjectOnNodalBasisUpdater<4>;
  template class ProjectOnNodalBasisUpdater<5>;
}
