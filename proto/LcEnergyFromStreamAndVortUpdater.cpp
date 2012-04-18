/**
 * @file	LcEnergyFromStreamAndVortUpdater.h
 *
 * @brief	Updater to compute energy from streamfunction.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcEnergyFromStreamAndVortUpdater.h>
#include <LcGlobals.h>
#include <LcLinAlgebra.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <vector>

namespace Lucee
{
  const char *EnergyFromStreamAndVortUpdater::id = "EnergyFromStreamAndVort";

  EnergyFromStreamAndVortUpdater::EnergyFromStreamAndVortUpdater()
    : Lucee::UpdaterIfc()
  {
  }
  
  void
  EnergyFromStreamAndVortUpdater::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::UpdaterIfc::readInput(tbl);

// get hold of element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis");
    else
      throw Lucee::Except("NodalFiniteElementIfc::readInput: Must specify element to use using 'basis'");
  }

  void
  EnergyFromStreamAndVortUpdater::initialize()
  {
// call base class method
    Lucee::UpdaterIfc::initialize();

// get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();

// local region to update
    Lucee::Region<2, int> localRgn = grid.getLocalRegion();

// set index to first location in grid (this is okay as in this
// updater we are assuming grid is uniform)
    nodalBasis->setIndex(localRgn.getLower(0), localRgn.getLower(1));

    unsigned nlocal = nodalBasis->getNumNodes();

// space for mass matrix
    Lucee::Matrix<double> massMatrix(nlocal, nlocal);

    for (unsigned dir=0; dir<2; ++dir)
    {
// get stiffness matrice
      Lucee::Matrix<double> stiffMatrix(nlocal, nlocal);
      nodalBasis->getGradStiffnessMatrix(dir, stiffMatrix);

// calculate differentiation matrix
      diffMatrix[dir].m = Lucee::Matrix<double>(nlocal, nlocal);
      for (unsigned i=0; i<nlocal; ++i)
        for (unsigned j=0; j<nlocal; ++j)
// diff matrices are computed from transposed stiffness matrices
          diffMatrix[dir].m(i,j) = stiffMatrix(j,i);

// multiply matrices by inverse of mass matrix
      nodalBasis->getMassMatrix(massMatrix);
      Lucee::solve(massMatrix, diffMatrix[dir].m);
    }
  }

  Lucee::UpdaterStatus
  EnergyFromStreamAndVortUpdater::update(double t)
  {
// get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();
// get input arrays
    const Lucee::Field<2, double>& phi = this->getInp<Lucee::Field<2, double> >(0);
    const Lucee::Field<2, double>& chi = this->getInp<Lucee::Field<2, double> >(1);
// get output dynVector
    Lucee::DynVector<double>& energy = this->getOut<Lucee::DynVector<double> >(0);

// local region to update
    Lucee::Region<2, int> localRgn = grid.getLocalRegion();
// global region
    Lucee::Region<2, int> globalRgn = grid.getGlobalRegion();

// number of local nodes
    unsigned nlocal = nodalBasis->getNumNodes();

// space for various quantities
    std::vector<double> phiK(nlocal), normGradPhi(nlocal);
    std::vector<double> weights(nlocal);
    Lucee::ConstFieldPtr<double> chiPtr = chi.createConstPtr();

    double totalEnergy = 0.0;
// compute total energy contribution from chi*phi term
    for (int ix=localRgn.getLower(0); ix<localRgn.getUpper(0); ++ix)
    {
      for (int iy=localRgn.getLower(1); iy<localRgn.getUpper(1); ++iy)
      {
        nodalBasis->setIndex(ix, iy);
        chi.setPtr(chiPtr, ix, iy);

// extract potential at this location
        nodalBasis->extractFromField(phi, phiK);

// get quadrature weights
        nodalBasis->getWeights(weights);

// compute contribution to energy from this cell
        for (unsigned k=0; k<nlocal; ++k)
          totalEnergy += weights[k]*chiPtr[k]*phiK[k];
      }
    }

    int idx[2];

// compute contribution from surface integrals to take into account
// "leakage" terms on UPPER surface
    for (unsigned d=0; d<2; ++d)
    {
// create region to loop over side
      Lucee::Region<2, int> defRgnG =
        globalRgn.resetBounds(d, globalRgn.getUpper(d)-1, globalRgn.getUpper(d)) ;
// only update if we are on the correct ranks
      Lucee::Region<2, int> defRgn = defRgnG.intersect(localRgn);

      unsigned nsurf = nodalBasis->getNumSurfUpperNodes(d);
// space for various things
      std::vector<double> grad(nlocal);
      std::vector<double> surfWeights(nsurf);
      std::vector<int> surfNodeNums(nsurf);

// get local node numbers on upper surface (get this once as it does
// not change)
      nodalBasis->getSurfUpperNodeNums(d, surfNodeNums);

      Lucee::RowMajorSequencer<2> seq(defRgn);
      while (seq.step())
      {
        seq.fillWithIndex(idx);
        nodalBasis->setIndex(idx);

// extract potential at this location
        nodalBasis->extractFromField(phi, phiK);
// compute gradient
        matVec(1.0, diffMatrix[d].m, phiK, 0.0, &grad[0]);

// get weights
        nodalBasis->getSurfUpperWeights(d, surfWeights);

// compute contribution to energy
        for (unsigned k=0; k<nsurf; ++k)
          totalEnergy += surfWeights[k]*phiK[surfNodeNums[k]-1]*grad[surfNodeNums[k]-1];
      }
    }

// compute contribution from surface integrals to take into account
// "leakage" terms on LOWER surface
    for (unsigned d=0; d<2; ++d)
    {
// create region to loop over side
      Lucee::Region<2, int> defRgnG =
        globalRgn.resetBounds(d, globalRgn.getLower(d), globalRgn.getLower(d)+1);
// only update if we are on the correct ranks
      Lucee::Region<2, int> defRgn = defRgnG.intersect(localRgn);

      unsigned nsurf = nodalBasis->getNumSurfLowerNodes(d);
// space for various things
      std::vector<double> grad(nlocal);
      std::vector<double> surfWeights(nsurf);
      std::vector<int> surfNodeNums(nsurf);

// get local node numbers on upper surface (get this once as it does
// not change)
      nodalBasis->getSurfLowerNodeNums(d, surfNodeNums);

      Lucee::RowMajorSequencer<2> seq(defRgn);
      while (seq.step())
      {
        seq.fillWithIndex(idx);
        nodalBasis->setIndex(idx);

// extract potential at this location
        nodalBasis->extractFromField(phi, phiK);
// compute gradient
        matVec(1.0, diffMatrix[d].m, phiK, 0.0, &grad[0]);

// get weights
        nodalBasis->getSurfLowerWeights(d, surfWeights);

// compute contribution to energy (lower surfaces contribute -ve as
// normal points inwards)
        for (unsigned k=0; k<nsurf; ++k)
          totalEnergy += -surfWeights[k]*phiK[surfNodeNums[k]-1]*grad[surfNodeNums[k]-1];
      }
    }

    double netTotalEnergy = totalEnergy;
// get hold of comm pointer to do all parallel messaging
    TxCommBase *comm = Loki::SingletonHolder<Lucee::Globals>
      ::Instance().comm;
// sum across all processors
    comm->allreduce(1, &totalEnergy, &netTotalEnergy, TX_SUM);

    std::vector<double> data(1);
    data[0] = 0.5*netTotalEnergy;
// push value into dynVector
    energy.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  void
  EnergyFromStreamAndVortUpdater::calcNormGrad(std::vector<double>& phiK,
    std::vector<double>& normGradPhi)
  {
// compute gradient in X- and Y-directions
    std::vector<double> gradX(normGradPhi.size()), gradY(normGradPhi.size());
    matVec(1.0, diffMatrix[0].m, phiK, 0.0, &gradX[0]);
    matVec(1.0, diffMatrix[1].m, phiK, 0.0, &gradY[0]);

// compute norm of gradient
    for (unsigned i=0; i<normGradPhi.size(); ++i)
      normGradPhi[i] = gradX[i]*gradX[i] + gradY[i]*gradY[i];
  }

  void
  EnergyFromStreamAndVortUpdater::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }

  void 
  EnergyFromStreamAndVortUpdater::matVec(double m, const Lucee::Matrix<double>& mat,
    const std::vector<double>& vec, double v, double *out)
  {
    double tv;
    unsigned rows = mat.numRows(), cols = mat.numColumns();
    for (unsigned i=0; i<rows; ++i)
    {
      tv = 0.0;
      for (unsigned j=0; j<cols; ++j)
        tv += mat(i,j)*vec[j];
      out[i] = m*tv + v*out[i];
    }
  }
}
