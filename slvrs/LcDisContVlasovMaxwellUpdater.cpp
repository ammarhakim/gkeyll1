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

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("DisContVlasovMaxwellUpdater::readInput: Must specify element to use using 'basis'");

    cfl = tbl.getNumber("cfl");
    cflm = 1.1*cfl; // use slightly large max CFL to avoid thrashing around

  }

  template <unsigned CDIM, unsigned VDIM>
  void
  DisContVlasovMaxwellUpdater<CDIM, VDIM>::initialize()
  {
    // call base class method
    UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    // local region to update
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    Lucee::RowMajorSequencer<CDIM> seq(localRgn);
    seq.step(); // just to get to first index
    int cidx[CDIM];
    seq.fillWithIndex(cidx);
    nodalBasis->setIndex(cidx);

    Lucee::RowMajorSequencer<VDIM> seq(localRgn);
    seq.step(); // just to get to first index
    int vidx[VDIM];
    seq.fillWithIndex(vidx);
    nodalBasis->setIndex(vidx);
    
    unsigned nlocal = nodalBasis->getNumNodes();

// get node numbers on each lower and upper edges
    for (unsigned dir=0; dir<NDIM; ++dir)
    {
      lowerNodeNums[dir].nums.resize(nodalBasis->getNumSurfLowerNodes(dir));
      nodalBasis->getSurfLowerNodeNums(dir, lowerNodeNums[dir].nums);

      upperNodeNums[dir].nums.resize(nodalBasis->getNumSurfUpperNodes(dir));
      nodalBasis->getSurfUpperNodeNums(dir, upperNodeNums[dir].nums);
    }

    Lucee::Matrix<double> massMatrix(nlocal, nlocal);

    for (unsigned dir=0; dir<NDIM; ++dir)
    {
      // get stiffness matrix
      stiffMatrix[dir].m = Lucee::Matrix<double>(nlocal, nlocal);
      nodalBasis->getGradStiffnessMatrix(dir, stiffMatrix[dir].m);

      nodalBasis->getMassMatrix(massMatrix);
      Lucee::solve(massMatrix, stiffMatrix[dir].m); // pre-multiply by inverse mass matrix

      // compute lift matrices
      lowerLift[dir].m = Lucee::Matrix<double>(nlocal, 
        nodalBasis->getNumSurfLowerNodes(dir));

      nodalBasis->getLowerFaceMassMatrix(dir, lowerLift[dir].m);
      nodalBasis->getMassMatrix(massMatrix);
      Lucee::solve(massMatrix, lowerLift[dir].m);  // pre-multiply by inverse mass matrix

      upperLift[dir].m = Lucee::Matrix<double>(nlocal, 
        nodalBasis->getNumSurfUpperNodes(dir));

      nodalBasis->getUpperFaceMassMatrix(dir, upperLift[dir].m);
      nodalBasis->getMassMatrix(massMatrix);
      Lucee::solve(massMatrix, upperLift[dir].m);  // pre-multiply by inverse mass matrix
    }

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
