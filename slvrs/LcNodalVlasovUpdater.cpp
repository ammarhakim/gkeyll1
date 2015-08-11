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
  bool
  NodalVlasovUpdater<CDIM,VDIM>::sameConfigCoords(unsigned n, unsigned cn, double dxMin,
    const Lucee::Matrix<double>& phaseC, const Lucee::Matrix<double>& confC)
  {
    for (unsigned d=0; d<CDIM; ++d)
      if (! (std::fabs(phaseC(n,d)-confC(cn,d))<1e-4*dxMin) )
        return false;
    return true;
  }

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
    confBasis->setIndex(idx); // only first CDIM elements are used
    
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

// compute mapping of phase-space nodes to configuration space
// nodes. The assumption here is that the node layout in phase-space
// and configuration space are such that each node in phase-space has
// exactly one node co-located with it in configuration space. No
// "orphan" phase-space node are allowed, and an exception is thrown
// if that occurs.
    phaseConfMap.resize(nlocal);
    Lucee::Matrix<double> phaseNodeCoords(phaseBasis->getNumNodes(), PNC);
    Lucee::Matrix<double> confNodeCoords(confBasis->getNumNodes(), CNC);

    double dxMin = grid.getDx(0);
    for (unsigned d=1; d<CDIM; ++d)
      dxMin = std::min(dxMin, grid.getDx(d));

    phaseBasis->getNodalCoordinates(phaseNodeCoords);
    confBasis->getNodalCoordinates(confNodeCoords);
    for (unsigned n=0; n<nlocal; ++n)
    {
      bool pcFound = false;
      for (unsigned cn=0; cn<nlocal; ++cn)
        if (sameConfigCoords(n, cn, dxMin, phaseNodeCoords, confNodeCoords))
        {
          phaseConfMap[n] = cn;
          pcFound = true;
          break;
        }
      if (!pcFound)
      {
        Lucee::Except lce(
          "NodalVlasovUpdater::readInput: No matching configuration space node for phase-space node ");
        lce << n;
        throw lce;
      }
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

      for (unsigned dir=0; dir<NDIM; ++dir)
      {
        for (unsigned n=0; n<nlocal; ++n)
        {
        }
      }
    }

// contributions from surface integrals
    for (unsigned dir=0; dir<NDIM; ++dir)
    {

// lower and upper bounds of 1D slice. (We need to make sure that flux
// is computed for one edge outside domain interior)
      int sliceLower = localRgn.getLower(dir);
      int sliceUpper = localRgn.getUpper(dir)+1;

      int idx[NDIM], idxl[NDIM];
// loop over each 1D slice
      while (seq.step())
      {
        seq.fillWithIndex(idx);
        seq.fillWithIndex(idxl);

        for (int i=sliceLower; i<sliceUpper; ++i)
        { // loop over each edge
          idx[dir] = i; // cell right of edge
          idxl[dir] = i-1; // cell left of edge

          grid.setIndex(idxl);
          double dxL = grid.getDx(dir);
          grid.setIndex(idx);
          double dxR = grid.getDx(dir);

          double dtdx = 2*dt/(dxL+dxR);

          q.setPtr(qPtr, idx);
          q.setPtr(qPtrl, idxl);

// compute numerical fluxes on each face
          unsigned nface = lowerNodeNums[dir].nums.size();

          for (unsigned s=0; s<nface; ++s)
          {
            unsigned un = upperNodeNums[dir].nums[s];
            unsigned ln = lowerNodeNums[dir].nums[s];

            //equation->rotateToLocal(coordSys, &qPtrl[meqn*un], &localQl[0]);
            //equation->rotateToLocal(coordSys, &qPtr[meqn*ln], &localQ[0]);

            double maxs = 0; //equation->numericalFlux(coordSys,
            //   &localQl[0], &localQ[0], inAuxQl, inAuxQr, &localF[0]);

            //equation->rotateToGlobal(coordSys, &localF[0], &flux[meqn*s]);

// compute actual CFL number to control time-stepping
            cfla = Lucee::max3(cfla, dtdx*maxs, -dtdx*maxs);
          }

// update left cell connected to edge with flux on face
          qNew.setPtr(qNewPtrl, idxl);
          matVec(-1.0, upperLift[dir].m, 1, &flux[0], 1.0, &qNewPtrl[0]);

// update right cell connected to edge with flux on face
          qNew.setPtr(qNewPtr, idx);
          matVec(1.0, lowerLift[dir].m, 1, &flux[0], 1.0, &qNewPtr[0]);
        }
      }
      if (cfla > cflm)
// time-step was too large: return a suggestion with correct time-step
        return Lucee::UpdaterStatus(false, dt*cfl/cfla);
    }

// NOTE: If only calculation of increments are requested, the final
// Euler update is not performed. This means that the multiplication
// of the DG RHS with dt is not done, something to keep in mind if
// using the increment in time-dependent update.
    if (onlyIncrement == false)
    {
      seq = Lucee::RowMajorSequencer<NDIM>(localRgn);
// final sweep, update solution with forward Euler step
      while (seq.step())
      {
        seq.fillWithIndex(idx);
        qNew.setPtr(qNewPtr, idx);
        q.setPtr(qPtr, idx);
        for (unsigned k=0; k<nlocal; ++k)
          qNewPtr[k] = qPtr[k] + dt*qNewPtr[k];
      }
    }

    return Lucee::UpdaterStatus(true, dt*cfl/cfla);
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
