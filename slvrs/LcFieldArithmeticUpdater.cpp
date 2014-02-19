/**
 * @file	LcFieldArithmeticUpdater.cpp
 *
 * @brief	Evaluate a function using provided (two) fields as input.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcFieldArithmeticUpdater.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <vector>

namespace Lucee
{

  template <> const char *FieldArithmeticUpdater<1>::id = "FieldArithmeticUpdater1D";
  template <> const char *FieldArithmeticUpdater<2>::id = "FieldArithmeticUpdater2D";
  template <> const char *FieldArithmeticUpdater<3>::id = "FieldArithmeticUpdater3D";

  template <unsigned NDIM>
  FieldArithmeticUpdater<NDIM>::FieldArithmeticUpdater()
    : fnRef(-1), sharedNodes(false)
  {
  }

  template <unsigned NDIM>
  void
  FieldArithmeticUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
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
      throw Lucee::Except("FieldArithmeticUpdater::readInput: Must specify element to use using 'basis'");
  }

  template <unsigned NDIM>
  void
  FieldArithmeticUpdater<NDIM>::initialize()
  {
    // call base class method
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  FieldArithmeticUpdater<NDIM>::update(double t)
  {
    // get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    // get input fields
    const Lucee::Field<NDIM, double>& field1In = this->getInp<Lucee::Field<NDIM, double> >(0);
    const Lucee::Field<NDIM, double>& field2In = this->getInp<Lucee::Field<NDIM, double> >(1);
    // get output array
    Lucee::Field<NDIM, double>& fieldOut = this->getOut<Lucee::Field<NDIM, double> >(0);
    
    fieldOut = 0.0; // set all entries to 0.0

    // Pointers
    Lucee::ConstFieldPtr<double> field1InPtr = field1In.createConstPtr();
    Lucee::ConstFieldPtr<double> field2InPtr = field2In.createConstPtr();
    Lucee::FieldPtr<double> fieldOutPtr = fieldOut.createPtr();

    // get list of nodes exclusively owned by element
    std::vector<int> ndIds;

    if (sharedNodes == true)
      nodalBasis->getExclusiveNodeIndices(ndIds);
    else
    {
      // create "unit" mapping
      ndIds.resize(nodalBasis->getNumNodes());

      for (unsigned i = 0; i<ndIds.size(); ++i)
        ndIds[i] = i;
    }

    // number of nodes
    unsigned numNodes = ndIds.size();
    // Needed variables in loop:
    int numLocalNodes = nodalBasis->getNumNodes();

    // determine number of components in field
    int nc = fieldOut.getNumComponents()/numNodes;
    // to store function evaluation result
    std::vector<double> resultVector(1); 
    // indices into grid
    int idx[NDIM];

    // get hold of Lua state object
    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;

    // loop over each cell in extended region
    Lucee::Region<NDIM, int> localExtRgn = fieldOut.getExtRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localExtRgn);

    while (seq.step())
    {
      seq.fillWithIndex(idx);
      field1In.setPtr(field1InPtr, idx);
      field2In.setPtr(field2InPtr, idx);
      fieldOut.setPtr(fieldOutPtr, idx);

      nodalBasis->setIndex(idx);
      
      // Loop over each node (and each component), calling the lua function
      // and then setting the value returned into output field
      for (int nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
      {
        for (int component = 0; component < nc; component++)
        {
          evaluateFunction(*L, t, field1InPtr[nc*component+nodeIndex],
            field2InPtr[nc*component+nodeIndex], resultVector);
          fieldOutPtr[nc*component+nodeIndex] = resultVector[0];
        }
      }
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  FieldArithmeticUpdater<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    // Output: some operation applied to the two input fields
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void
  FieldArithmeticUpdater<NDIM>::evaluateFunction(Lucee::LuaState& L, double tm,
    double field1Component, double field2Component, std::vector<double>& res)
  {
    // push function object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, fnRef);
    // push variables on stack
    lua_pushnumber(L, field1Component);
    lua_pushnumber(L, field2Component);
    lua_pushnumber(L, tm);
    // call function
    if (lua_pcall(L, 3, res.size(), 0) != 0)
    {
      Lucee::Except lce("FieldArithmeticUpdater::evaluateFunction: ");
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
        throw Lucee::Except("FieldArithmeticUpdater::evaluateFunction: Return value not a number");
      res[res.size()+i] = lua_tonumber(L, i);
    }
    lua_pop(L, 1);
  }

  // instantiations
  template class FieldArithmeticUpdater<1>;
  template class FieldArithmeticUpdater<2>;
  template class FieldArithmeticUpdater<3>;
}
