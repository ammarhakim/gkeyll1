/**
 * @file	 DisContVlasovMaxwellUpdater.cpp
 *
 * @brief        Vlasov-Maxwell solver
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDisContVlasovMaxwellUpdater.h>

namespace Lucee
{

// set ids for module system
  template <> const char* DisContVlasovMaxwellUpdater<1,1>::id = "DgVlasovMaxwell1X1V";
  template <> const char* DisContVlasovMaxwellUpdater<1,2>::id = "DgVlasovMaxwell1X2V";
  template <> const char* DisContVlasovMaxwellUpdater<1,3>::id = "DgVlasovMaxwell1X3V";
  template <> const char* DisContVlasovMaxwellUpdater<2,2>::id = "DgVlasovMaxwell2X2V";
  template <> const char* DisContVlasovMaxwellUpdater<2,3>::id = "DgVlasovMaxwell2X3V";
  template <> const char* DisContVlasovMaxwellUpdater<3,3>::id = "DgVlasovMaxwell3X3V";

  template <unsigned CDIM, unsigned VDIM>
  DisContVlasovMaxwellUpdater<CDIM, VDIM>::~DisContVlasovMaxwellUpdater()
  {

  }

  template <unsigned CDIM, unsigned VDIM>
  void
  DisContVlasovMaxwellUpdater<CDIM, VDIM>::declareTypes()
  {

  }

  template <unsigned CDIM, unsigned VDIM>
  void
  DisContVlasovMaxwellUpdater<CDIM, VDIM>::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    UpdaterIfc::readInput(tbl);
  }

  template <unsigned CDIM, unsigned VDIM>
  void
  DisContVlasovMaxwellUpdater<CDIM, VDIM>::initialize()
  {
    // call base class method
    UpdaterIfc::initialize();
  }

  template <unsigned CDIM, unsigned VDIM> 
  Lucee::UpdaterStatus 
  DisContVlasovMaxwellUpdater<CDIM, VDIM>::update(double t)
  {

     return UpdaterStatus();
  }

// instantiations
  template class DisContVlasovMaxwellUpdater<1,1>;
  template class DisContVlasovMaxwellUpdater<1,2>;
  template class DisContVlasovMaxwellUpdater<1,3>;
  template class DisContVlasovMaxwellUpdater<2,2>;
  template class DisContVlasovMaxwellUpdater<2,3>;
  template class DisContVlasovMaxwellUpdater<3,3>;
}
