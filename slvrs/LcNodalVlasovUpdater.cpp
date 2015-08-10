/**
 * @file	LcNodalVlasovUpdater.cpp
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
#include <LcNodalVlasovUpdater.h>
#include <LcStructuredGridBase.h>

// std includes
#include <vector>

namespace Lucee
{
// set id for module system
  template <> const char *NodalVlasovUpdater<1,1>::id = "NodalVlasov1X1V";
  template <> const char *NodalVlasovUpdater<1,2>::id = "NodalVlasov1X2V";
  template <> const char *NodalVlasovUpdater<1,3>::id = "NodalVlasov1X3V";
  template <> const char *NodalVlasovUpdater<2,2>::id = "NodalVlasov2X2V";
  template <> const char *NodalVlasovUpdater<2,3>::id = "NodalVlasov2X3V";
  //template <> const char *NodalVlasovUpdater<3,3>::id = "NodalVlasov3X3V";

  template <unsigned CDIM, unsigned VDIM>
  NodalVlasovUpdater<CDIM,VDIM>::NodalVlasovUpdater()
    : UpdaterIfc()
  {
  }

  template <unsigned CDIM, unsigned VDIM>  
  void 
  NodalVlasovUpdater<CDIM,VDIM>::readInput(Lucee::LuaTable& tbl)
  {
    const unsigned NDIM = CDIM+VDIM;
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("phaseBasis"))
      phaseBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("phaseBasis");
    else
      throw Lucee::Except("NodalVlasovUpdater::readInput: Must specify phase-space basis using 'phaseBasis'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<CDIM> >("confBasis"))
      confBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<CDIM> >("confBasis");
    else
      throw Lucee::Except("NodalVlasovUpdater::readInput: Must specify configuration-space basis using 'confBasis'");

// directions to update
    if (tbl.hasNumVec("updateDirections"))
    {
      std::vector<double> ud = tbl.getNumVec("updateDirections");
      for (unsigned i=0; i<ud.size(); ++i)
      {
        unsigned d = (unsigned) ud[i];
        if (d<NDIM)
          updateDims.push_back(d);
        else
        {
          Lucee::Except lce("updateDirections must be a table less than ");
          lce << NDIM;
          throw lce;
        }
      }
    }
    else
    {
      for (unsigned i=0; i<NDIM; ++i)
        updateDims.push_back(i);
    }

    cfl = tbl.getNumber("cfl");
    cflm = 1.1*cfl; // use slightly large max CFL to avoid thrashing around

    onlyIncrement = false;
// when onlyIncrement flag is set contribution is not added to the
// input field, i.e. only increment is computed
    if (tbl.hasBool("onlyIncrement"))
      onlyIncrement = tbl.getBool("onlyIncrement");
  }

  template <unsigned CDIM, unsigned VDIM>
  void 
  NodalVlasovUpdater<CDIM,VDIM>::initialize()
  {
    const unsigned NDIM = CDIM+VDIM;
    Lucee::UpdaterIfc::initialize();

// get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
// local region to update
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[NDIM];
    seq.fillWithIndex(idx);
    phaseBasis->setIndex(idx);
    
    unsigned nlocal = phaseBasis->getNumNodes();

// get node numbers on each lower and upper edges
    for (unsigned dir=0; dir<NDIM; ++dir)
    {
      lowerNodeNums[dir].nums.resize(phaseBasis->getNumSurfLowerNodes(dir));
      phaseBasis->getSurfLowerNodeNums(dir, lowerNodeNums[dir].nums);

      upperNodeNums[dir].nums.resize(phaseBasis->getNumSurfUpperNodes(dir));
      phaseBasis->getSurfUpperNodeNums(dir, upperNodeNums[dir].nums);
    }

    Lucee::Matrix<double> massMatrix(nlocal, nlocal);

    for (unsigned dir=0; dir<NDIM; ++dir)
    {
// get stiffness matrix
      stiffMatrix[dir].m = Lucee::Matrix<double>(nlocal, nlocal);
      phaseBasis->getGradStiffnessMatrix(dir, stiffMatrix[dir].m);

      phaseBasis->getMassMatrix(massMatrix);
      Lucee::solve(massMatrix, stiffMatrix[dir].m); // pre-multiply by inverse mass matrix

// compute lift matrices
      lowerLift[dir].m = Lucee::Matrix<double>(nlocal, 
        phaseBasis->getNumSurfLowerNodes(dir));

      phaseBasis->getLowerFaceMassMatrix(dir, lowerLift[dir].m);
      phaseBasis->getMassMatrix(massMatrix);
      Lucee::solve(massMatrix, lowerLift[dir].m);  // pre-multiply by inverse mass matrix

      upperLift[dir].m = Lucee::Matrix<double>(nlocal, 
        phaseBasis->getNumSurfUpperNodes(dir));

      phaseBasis->getUpperFaceMassMatrix(dir, upperLift[dir].m);
      phaseBasis->getMassMatrix(massMatrix);
      Lucee::solve(massMatrix, upperLift[dir].m);  // pre-multiply by inverse mass matrix
    }    
  }

  template <unsigned CDIM, unsigned VDIM>
  Lucee::UpdaterStatus 
  NodalVlasovUpdater<CDIM,VDIM>::update(double t)
  {
    const unsigned NDIM = CDIM+VDIM;    
    const Lucee::StructuredGridBase<NDIM>& grid
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

// for compatibility with NodalDisContHyperUpdater I am retaining the
// names "q" and "qNew" for the distribution function. Hence, q is the
// distribution function at time t and qNew the distribution function
// at t+dt (Ammar Hakim)
    const Lucee::Field<NDIM, double>& q = this->getInp<Lucee::Field<NDIM, double> >(0);
    const Lucee::Field<CDIM, double>& EM = this->getInp<Lucee::Field<CDIM, double> >(1);
    Lucee::Field<NDIM, double>& qNew = this->getOut<Lucee::Field<NDIM, double> >(0);

    unsigned nlocal = phaseBasis->getNumNodes();

    double dt = t-this->getCurrTime();
    double cfla = 0.0; // maximum CFL number used    
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    Lucee::ConstFieldPtr<double> qPtr = q.createConstPtr();
    Lucee::ConstFieldPtr<double> qPtrl = q.createConstPtr();
    Lucee::FieldPtr<double> qNewPtr = qNew.createPtr();
    Lucee::FieldPtr<double> qNewPtrl = qNew.createPtr();
    std::vector<double> flux(nlocal);
    double localQ, localQl, localF;

    qNew = 0.0; // use qNew to store increment initially    
    
    int idx[NDIM];
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);

// loop to compute contribution from volume integrals
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      q.setPtr(qPtr, idx);
      qNew.setPtr(qNewPtr, idx);

      for (unsigned d=0; d<updateDims.size(); ++d)
      {
        unsigned dir = updateDims[d]; // direction to update
        for (unsigned n=0; n<nlocal; ++n)
        {
        }
        //matVec(1.0, stiffMatrix[dir].m, 1, &flux[0], 1.0, &qNewPtr[0]); // stiffness X flux
      }          
    }
    
    return Lucee::UpdaterStatus();
  }

  template <unsigned CDIM, unsigned VDIM>  
  void
  NodalVlasovUpdater<CDIM,VDIM>::declareTypes()
  {
    const unsigned NDIM = CDIM+VDIM;    
// distribution function
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
// E and B field in a single field
    this->appendInpVarType(typeid(Lucee::Field<CDIM, double>));
// returns one output: updated distribution function
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned CDIM, unsigned VDIM>
  void 
  NodalVlasovUpdater<CDIM,VDIM>::matVec(double mc, const Lucee::Matrix<double>& mat,
    unsigned meqn, const double* vec, double v, double *out)
  {
    double tv;
    unsigned rows = mat.numRows(), cols = mat.numColumns();
    for (unsigned m=0; m<meqn; ++m)
    {
      for (unsigned i=0; i<rows; ++i)
      {
        tv = 0.0;
        for (unsigned j=0; j<cols; ++j)
          tv += mat(i,j)*vec[meqn*j+m];
        out[meqn*i+m] = mc*tv + v*out[meqn*i+m];
      }
    }
  }

// instantiations
  template class NodalVlasovUpdater<1,1>;
  template class NodalVlasovUpdater<1,2>;
  template class NodalVlasovUpdater<1,3>;
  template class NodalVlasovUpdater<2,2>;
  template class NodalVlasovUpdater<2,3>;
  //template class NodalVlasovUpdater<3,3>;
}
