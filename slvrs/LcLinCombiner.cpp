/**
 * @file	LcLinCombiner.cpp
 *
 * @brief	Linear combiner updater.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcLinCombiner.h>
#include <LcStructuredGridBase.h>

namespace Lucee
{
// set ids for creators
  template <> const char *LinCombiner<1>::id = "LinCombiner1D";
  template <> const char *LinCombiner<2>::id = "LinCombiner2D";
  template <> const char *LinCombiner<3>::id = "LinCombiner3D";

  template <unsigned NDIM>
  LinCombiner<NDIM>::LinCombiner()
    : Lucee::UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  LinCombiner<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    UpdaterIfc::readInput(tbl);

    if (tbl.hasNumVec("coeffs"))
      coeff = tbl.getNumVec("coeffs");
    else
      throw Lucee::Except("LinCombiner<NDIM>::readInput: Must provide a 'coeffs' array.");
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  LinCombiner<NDIM>::update(double t)
  {
// check if number of coeffs matched number of input data structures
    if (coeff.size() != this->getNumInpVars())
    {
      Lucee::Except lce("LinCombiner<NDIM>::update: Number of coefficients (");
      lce << coeff.size() << ") should be same as number of input fields ("
          << this->getNumInpVars() << ")";
      throw lce;
    }

// get output field
    Lucee::Field<NDIM, double>& outFld = this->template 
      getOut<Lucee::Field<NDIM, double> >(0);

// set output array to 0.0
    outFld = 0.0;

// combine input arrays accumulating it in output array
    for (unsigned i=0; i<this->getNumInpVars(); ++i)
    {
      const Lucee::Field<NDIM, double>& inpFld = this->template 
        getInp<Lucee::Field<NDIM, double> >(i);
      outFld.accumulate(coeff[i], inpFld);
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  LinCombiner<NDIM>::declareTypes()
  {
// only one output field is expected
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
// arbitrary number of input fields
    this->setLastInpVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void
  LinCombiner<NDIM>::setCoeff(const std::vector<double>& c)
  {
    coeff = c;
  }

// instantiations
  template class LinCombiner<1>;
  template class LinCombiner<2>;
  template class LinCombiner<3>;
}
