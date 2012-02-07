/**
 * @file	LcProjectOnBasisUpdater.cpp
 *
 * @brief	Project a function of a basis functions.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGlobals.h>
#include <LcMathLib.h>
#include <LcProjectOnBasisUpdater.h>
#include <LcStructuredGridBase.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  template <> const char *ProjectOnBasisUpdater<1>::id = "ProjectOnBasis1D";
  template <> const char *ProjectOnBasisUpdater<2>::id = "ProjectOnBasis2D";
  template <> const char *ProjectOnBasisUpdater<3>::id = "ProjectOnBasis3D";

  template <unsigned NDIM>
  ProjectOnBasisUpdater<NDIM>::ProjectOnBasisUpdater()
    : numBasis(1), fnRef(-1), Pmk(1,1), w(1), mu(1)
  {
// this ctor should not be used directly
  }

  template <unsigned NDIM>
  void
  ProjectOnBasisUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    UpdaterIfc::readInput(tbl);
// get number of basis functions to project on
    numBasis = (unsigned) tbl.getNumber("numBasis");
// get function to project
    fnRef = tbl.getFunctionRef("project");

// allocate space
    Pmk = Matrix<double>(numBasis, numBasis);
    w = Lucee::Vector<double>(numBasis);
    mu = Lucee::Vector<double>(numBasis);

// compute weights and ordinates
    Lucee::gauleg(numBasis, -1, 1, mu, w);

// compute Legendre polynomials at ordinates
    for (unsigned m=0; m<numBasis; ++m)
      for (unsigned k=0; k<numBasis; ++k)
        Pmk(m,k) = Lucee::legendrePoly(m, mu[k]);
  }

  template <unsigned NDIM>
  void
  ProjectOnBasisUpdater<NDIM>::initialize()
  {
// call base class method
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  ProjectOnBasisUpdater<NDIM>::update(double t)
  {
// get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
// get output array
    Lucee::Field<NDIM, double>& q = this->getOut<Lucee::Field<NDIM, double> >(0);

// create region to walk over ordinates
    int shape[NDIM];
    for (unsigned i=0; i<NDIM; ++i) shape[i] = numBasis;
    Lucee::Region<NDIM, int> ordRgn(shape);
// create sequencer and indexer on this region
    Lucee::RowMajorSequencer<NDIM> ordSeq(ordRgn);
    Lucee::RowMajorIndexer<NDIM> ordIdxr(ordRgn);

// determine number of components in field (THIS IS NOT GOOD PRACTICE
// AND NEEDS TO CHANGE EVENTUALLY. SOMEHOW NODAL LAYOUT NEEDS TO BE
// FIELD OR SOMEWHERE BETTER. PRESENTLY THE NODAL LAYOUT IS STORED
// ONLY IMPLICITLY BY CONTRACT. Ammar Hakim, Feb 7 2012).
    unsigned nc = q.getNumComponents()/ordRgn.getVolume();

// indices into grid and ordinate in cell
    int idx[NDIM], ordIdx[NDIM];
// coordinates of cell centroid and ordinate in cell
    double xc[3], xmu[3];

// get hold of Lua state object
    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;

// pointer to data in field
    Lucee::FieldPtr<double> ptr = q.createPtr();

// loop over each cell in extended region
    Lucee::Region<NDIM, int> localExtRgn = q.getExtRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localExtRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      q.setPtr(ptr, idx);

// compute cell center coordinates
      grid.setIndex(idx);
      grid.getCentroid(xc);


    }
  }

  template <unsigned NDIM>
  void
  ProjectOnBasisUpdater<NDIM>::declareTypes()
  {
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void
  ProjectOnBasisUpdater<NDIM>::evaluateFunction(
    Lucee::LuaState& L, double tm, const double loc[3], std::vector<double>& res)
  {
// push function object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, fnRef);
// push variables on stack
    for (unsigned i=0; i<3; ++i)
      lua_pushnumber(L, loc[i]);
    lua_pushnumber(L, tm);
// call function
    if (lua_pcall(L, 4, res.size(), 0) != 0)
    {
      Lucee::Except lce("ProjectOnBasisUpdater::evaluateFunction: ");
      lce << "Problem evaluating function supplied as 'project' ";
      throw lce;
    }
// fetch results
    for (int i=-res.size(); i<0; ++i)
    {
      if (!lua_isnumber(L, i))
        throw Lucee::Except("ProjectOnBasisUpdater::evaluateFunction::getSource: Return value not a number");
      res[res.size()+i] = lua_tonumber(L, i);
    }
    lua_pop(L, 1);
  }

// instantiations
  template class ProjectOnBasisUpdater<1>;
  template class ProjectOnBasisUpdater<2>;
  template class ProjectOnBasisUpdater<3>;
}
