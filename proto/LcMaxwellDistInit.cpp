/**
 * @file	LcMaxwellDistInit.cpp
 *
 * @brief	Updater for initializin the Maxwellian distribution from moments
 */

// lucee includes
#include <LcMaxwellDistInit.h>
#include <LcStructGridField.h>

// eigen inlcudes
#include <Eigen/Eigen>

namespace Lucee
{
// set ids for module system
  template <> const char *MaxwellDistInit<1>::id = "MaxwellDistInit1D";
  template <> const char *MaxwellDistInit<2>::id = "MaxwellDistInit2D";
  template <> const char *MaxwellDistInit<3>::id = "MaxwellDistInit3D";

  template <unsigned NDIM>
  void
  MaxwellDistInit<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);
    // dir = (unsigned) tbl.getNumber("dir");
    // gravity = tbl.getNumber("gravity");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("NodalDisContSrcIncrUpdater::readInput: Must specify element to use using 'basis'");
  }

  template <unsigned NDIM>
  void
  MaxwellDistInit<NDIM>::initialize()
  {
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  MaxwellDistInit<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();

    double dt = t-this->getCurrTime();
    //    Lucee::Field<2, double>& fluid = this->getOut<Lucee::Field<2, double> >(0);
    //    Lucee::FieldPtr<double> fPtr = fluid.createPtr();
    
    return Lucee::UpdaterStatus();
  }
  
  template <unsigned NDIM>
  void
  MaxwellDistInit<NDIM>::declareTypes()
  {
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }

// instantiations
  template class MaxwellDistInit<1>;
  template class MaxwellDistInit<2>;
  template class MaxwellDistInit<3>;
}

