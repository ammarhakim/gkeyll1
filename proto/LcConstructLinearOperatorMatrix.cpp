/**
 * @file	LcConstructLinearOperatorMatrix.cpp
 *
 * @brief	Set region based on predicate.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcConstructLinearOperatorMatrix.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <vector>
// eigen include to write to file
#include <unsupported/Eigen/SparseExtra>

namespace Lucee
{

  template <> const char *ConstructLinearOperatorMatrix<1>::id = "ConstructLinearOperatorMatrix1D";
  template <> const char *ConstructLinearOperatorMatrix<2>::id = "ConstructLinearOperatorMatrix2D";
  template <> const char *ConstructLinearOperatorMatrix<3>::id = "ConstructLinearOperatorMatrix3D";
  template <> const char *ConstructLinearOperatorMatrix<4>::id = "ConstructLinearOperatorMatrix4D";
  template <> const char *ConstructLinearOperatorMatrix<5>::id = "ConstructLinearOperatorMatrix5D";

  template <unsigned NDIM>
  ConstructLinearOperatorMatrix<NDIM>::ConstructLinearOperatorMatrix()
    : UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  ConstructLinearOperatorMatrix<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);
    // get function to evaluate
    fnRef = tbl.getFunctionRef("evaluate");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("ConstructLinearOperatorMatrix::readInput: Must specify element to use using 'basis'");

    if (tbl.hasNumber("totalNodes"))
      totalNodes = tbl.getNumber("totalNodes");
    else
      throw Lucee::Except("ConstructLinearOperatorMatrix::readInput: Must specify totalNodes");
  }

  template <unsigned NDIM>
  void
  ConstructLinearOperatorMatrix<NDIM>::initialize()
  {
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  ConstructLinearOperatorMatrix<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    // get input field
    const Lucee::Field<NDIM, double>& fldIn = this->getInp<Lucee::Field<NDIM, double> >(0);

    std::vector<double> res(1);
    int idx[NDIM];
    int nlocal = nodalBasis->getNumNodes();
    
    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;

    Lucee::ConstFieldPtr<double> fldPtr = fldIn.createConstPtr();
    
    Lucee::Region<NDIM, int> localExtRgn = fldIn.getExtRegion();
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localExtRgn);
    Lucee::RowMajorIndexer<NDIM> volIdxr(localExtRgn);

    // Figure out what column to write to
    evaluateFunction(*L, t, res);
    int targetColumn = (int) res[0];

    // Loop over entirety of fldIn, storing only non-zero elements
    while(seq.step())
    {
      seq.fillWithIndex(idx);
      fldIn.setPtr(fldPtr, idx);
      nodalBasis->setIndex(idx);

      int cellIndex = volIdxr.getIndex(idx);

      for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
      {
        if(fldPtr[nodeIndex] != 0.0)
          tripletList.push_back(Eigen::Triplet<double>
            (cellIndex*nlocal+nodeIndex,targetColumn,fldPtr[nodeIndex]));
      }
    }

    if (targetColumn == totalNodes-1)
    {
      // Create sparse matrix from list of triplets
      Eigen::SparseMatrix<double> linearOperatorMatrix(totalNodes, totalNodes);
      linearOperatorMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
      // Save matrix to disk
      Eigen::saveMarket(linearOperatorMatrix,"linearOperatorMatrix.mtx");
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  ConstructLinearOperatorMatrix<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void
  ConstructLinearOperatorMatrix<NDIM>::evaluateFunction(Lucee::LuaState& L, double tm,
    std::vector<double>& res)
  {
// push function object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, fnRef);
// push variables on stack
    lua_pushnumber(L, tm);
// call function
    if (lua_pcall(L, 1, res.size(), 0) != 0)
    {
      Lucee::Except lce("ConstructLinearOperatorMatrix::evaluateFunction: ");
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
        throw Lucee::Except("ConstructLinearOperatorMatrix::evaluateFunction: Return value not a number");
      res[res.size()+i] = lua_tonumber(L, i);
    }
    lua_pop(L, 1);
  }

// instantiations
  template class ConstructLinearOperatorMatrix<1>;
  template class ConstructLinearOperatorMatrix<2>;
  template class ConstructLinearOperatorMatrix<3>;
  template class ConstructLinearOperatorMatrix<4>;
  template class ConstructLinearOperatorMatrix<5>;
}
