/**
 * @file	LcNodalDisContHyperUpdater.cpp
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
#include <LcNodalDisContHyperUpdater.h>
#include <LcStructuredGridBase.h>

// std includes
#include <vector>

namespace Lucee
{
// set id for module system
  template <> const char *NodalDisContHyperUpdater<1>::id = "NodalDgHyper1D";
  template <> const char *NodalDisContHyperUpdater<2>::id = "NodalDgHyper2D";
  template <> const char *NodalDisContHyperUpdater<3>::id = "NodalDgHyper3D";

  template <unsigned NDIM>
  NodalDisContHyperUpdater<NDIM>::NodalDisContHyperUpdater()
    : UpdaterIfc()
  {
  }

  template <unsigned NDIM>  
  void 
  NodalDisContHyperUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("NodalDisContHyperUpdater::readInput: Must specify element to use using 'basis'");

    if (tbl.hasObject<Lucee::HyperEquation>("equation"))
      equation = &tbl.getObjectAsBase<Lucee::HyperEquation>("equation");
    else
    {
      Lucee::Except lce("NodalDisContHyperUpdater::readInput: Must specify an equation to solve!");
      throw lce;
    }

    cfl = tbl.getNumber("cfl");
    cflm = 1.1*cfl; // use slightly large max CFL to avoid thrashing around

// when onlyIncrement flag is set contribution is not added to the
// input field, i.e. only increment is computed
    onlyIncrement = false;
    if (tbl.hasBool("onlyIncrement"))
      onlyIncrement = tbl.getBool("onlyIncrement");
  }

  template <unsigned NDIM>
  void 
  NodalDisContHyperUpdater<NDIM>::initialize()
  {
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
    nodalBasis->setIndex(idx);
    
    unsigned nlocal = nodalBasis->getNumNodes();

// get node numbers on each lower and upper edges
    for (unsigned dir=0; dir<NDIM; ++dir)
    {
      lowerNodeNums[dir].nums.resize(nodalBasis->getNumSurfLowerNodes(dir));
      nodalBasis->getSurfLowerNodeNums(dir, lowerNodeNums[dir].nums);

// reset numbers as element offsets them with 1
      for (unsigned k=0; k<nodalBasis->getNumSurfLowerNodes(dir); ++k)
        lowerNodeNums[dir].nums[k] += -1;

      upperNodeNums[dir].nums.resize(nodalBasis->getNumSurfUpperNodes(dir));
      nodalBasis->getSurfUpperNodeNums(dir, upperNodeNums[dir].nums);

// reset numbers as element offsets them with 1
      for (unsigned k=0; k<nodalBasis->getNumSurfUpperNodes(dir); ++k)
        upperNodeNums[dir].nums[k] += -1;
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

  template <unsigned NDIM>
  Lucee::UpdaterStatus 
  NodalDisContHyperUpdater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    const Lucee::Field<NDIM, double>& q = this->getInp<Lucee::Field<NDIM, double> >(0);
    Lucee::Field<NDIM, double>& qNew = this->getOut<Lucee::Field<NDIM, double> >(0);

    unsigned nlocal = nodalBasis->getNumNodes();
    unsigned meqn = equation->getNumEqns();

// determine number of auxillary variables
    unsigned numAuxVars = this->getNumInpVars()-1; // first is alway conserved variable
// store them
    std::vector<const Lucee::Field<NDIM, double>* > auxVars;
    for (unsigned i=0; i<numAuxVars; ++i)
      auxVars.push_back( &this->getInp<Lucee::Field<NDIM, double> >(i+1) );
// compute number of equations in each auxillary variable
    std::vector<unsigned> numAuxEqns;
    for (unsigned i=0; i<numAuxVars; ++i)
      numAuxEqns.push_back( auxVars[i]->getNumComponents()/nlocal );

    std::vector<const double *> auxQ(numAuxVars), auxQl(numAuxVars), auxQr(numAuxVars);
    std::vector<const double *> inAuxQl(numAuxVars), inAuxQr(numAuxVars);

    double dt = t-this->getCurrTime();
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    double cfla = 0.0; // maximum CFL number used

    Lucee::ConstFieldPtr<double> qPtr = q.createConstPtr();
    Lucee::ConstFieldPtr<double> qPtrl = q.createConstPtr();
    Lucee::FieldPtr<double> qNewPtr = qNew.createPtr();
    Lucee::FieldPtr<double> qNewPtrl = qNew.createPtr();
    std::vector<double> flux(nlocal*meqn);
    std::vector<double> localQ(meqn), localQl(meqn), localF(meqn);

    qNew = 0.0; // use qNew to store increment initially

    int idx[NDIM];
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);

// loop to compute contribution from volume integrals
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      q.setPtr(qPtr, idx);
      qNew.setPtr(qNewPtr, idx);

// create input auxiliary variables list
      for (unsigned a=0; a<numAuxVars; ++a)
      {
        Lucee::ConstFieldPtr<double> aPtr = auxVars[a]->createConstPtr();
        auxVars[a]->setPtr(aPtr, idx);
        auxQ[a] = &aPtr[0];
      }

      for (unsigned dir=0; dir<NDIM; ++dir)
      {
// volume fluxes should always use aligned CS, even on non-rectangular geometery
        Lucee::AlignedRectCoordSys coordSys(dir);
        for (unsigned n=0; n<nlocal; ++n)
        {
          equation->rotateToLocal(coordSys, &qPtr[meqn*n], &localQ[0]);
          equation->flux(coordSys, &localQ[0], auxQ, &localF[0]);
          equation->rotateToGlobal(coordSys, &localF[0], &flux[meqn*n]);
//  adjust auxilary variable pointers to point to data at next node
          for (unsigned a=0; a<numAuxVars; ++a)
            auxQ[a] = auxQ[a] + numAuxEqns[a];
        }
        matVec(1.0, stiffMatrix[dir].m, meqn, &flux[0], 1.0, &qNewPtr[0]); // stiffness X flux
      }
    }

// contributions from surface integrals
    for (unsigned dir=0; dir<NDIM; ++dir)
    {
      Lucee::AlignedRectCoordSys coordSys(dir); // eventually this needs to be a CS on the face
// create sequencer to loop over *each* 1D slice in 'dir' direction
      Lucee::RowMajorSequencer<NDIM> seq(localRgn.deflate(dir));

// lower and upper bounds of 1D slice. (We need to make sure that the
// flux is computed for one edge outside the domain interior)
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

// create input auxiliary variables list
          for (unsigned a=0; a<numAuxVars; ++a)
          {
            Lucee::ConstFieldPtr<double> aPtr = auxVars[a]->createConstPtr();
            auxVars[a]->setPtr(aPtr, idxl);
            auxQl[a] = &aPtr[0];

            auxVars[a]->setPtr(aPtr, idx);
            auxQr[a] = &aPtr[0];
          }

// compute numerical fluxes on each face
          unsigned nface = lowerNodeNums[dir].nums.size();
          for (unsigned s=0; s<nface; ++s)
          {
            unsigned un = upperNodeNums[dir].nums[s];
            unsigned ln = lowerNodeNums[dir].nums[s];

            equation->rotateToLocal(coordSys, &qPtrl[meqn*un], &localQl[0]);
            equation->rotateToLocal(coordSys, &qPtr[meqn*ln], &localQ[0]);
// NOTE: we do not rotate auxiliary variables as the equation system
// should do the rotations if needed. This is perhaps inconsistent,
// but the HyperEquation class interface need not be cluttered with
// yet another set of rotation functions.

// adjust auxilary variable pointers to point to data at face node
            for (unsigned a=0; a<numAuxVars; ++a)
            {
              inAuxQl[a] = auxQl[a] + un*numAuxEqns[a];
              inAuxQr[a] = auxQr[a] + ln*numAuxEqns[a];
            }

            double maxs = equation->numericalFlux(coordSys,
              &localQl[0], &localQ[0], inAuxQl, inAuxQr, &localF[0]);

            equation->rotateToGlobal(coordSys, &localF[0], &flux[meqn*s]);

// compute actual CFL number to control time-stepping
            cfla = Lucee::max3(cfla, dtdx*maxs, -dtdx*maxs);
          }

// update left cell connected to edge with flux on face
          qNew.setPtr(qNewPtrl, idxl);
          matVec(-1.0, upperLift[dir].m, meqn, &flux[0], 1.0, &qNewPtrl[0]);

// update right cell connected to edge with flux on face
          qNew.setPtr(qNewPtr, idx);
          matVec(1.0, lowerLift[dir].m, meqn, &flux[0], 1.0, &qNewPtr[0]);
        }
      }
      if (cfla > cflm)
// time-step was too large: return a suggestion with correct time-step
        return Lucee::UpdaterStatus(false, dt*cfl/cfla);
    }

    if (onlyIncrement == false)
    {
// NOTE: If only calculation of increments are requested, the final
// Euler update is not performed. This means that the multiplication
// of the DG RHS with dt is not done, something to keep in mind if
// using the increment in time-dependent update.

      seq = Lucee::RowMajorSequencer<NDIM>(localRgn);
// final sweep, update solution with forward Euler step
      while (seq.step())
      {
        seq.fillWithIndex(idx);
        qNew.setPtr(qNewPtr, idx);
        q.setPtr(qPtr, idx);
        
        for (unsigned k=0; k<qPtr.getNumComponents(); ++k)
          qNewPtr[k] = qPtr[k] + dt*qNewPtr[k];
      }
    }

    return Lucee::UpdaterStatus(true, dt*cfl/cfla);
  }

  template <unsigned NDIM>  
  void
  NodalDisContHyperUpdater<NDIM>::declareTypes()
  {
// takes one input
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
// returns one output
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void 
  NodalDisContHyperUpdater<NDIM>::matVec(double mc, const Lucee::Matrix<double>& mat,
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
  template class NodalDisContHyperUpdater<1>;
  template class NodalDisContHyperUpdater<2>;
  template class NodalDisContHyperUpdater<3>;
}
