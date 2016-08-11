/**
 * @file	LcNonLinEmGke1dHamilUpdater.cpp
 *
 * @brief	Compute nonlinear Hamiltonian for 1D EM/GKE
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcNonLinEmGke1dHamilUpdater.h>

namespace Lucee
{
  const char *NonLinEmGke1dHamilUpdater::id = "NonLinEmGke1DHamil";

  NonLinEmGke1dHamilUpdater::NonLinEmGke1dHamilUpdater()
    : UpdaterIfc()
  {
  }

  void
  NonLinEmGke1dHamilUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);
    charge = tbl.getNumber("charge");
    mass = tbl.getNumber("mass");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis");
    else
      throw Lucee::Except("NonLinEmGke1dHamilUpdater::readInput: Must specify element to use using 'basis'");
  }

  void
  NonLinEmGke1dHamilUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus
  NonLinEmGke1dHamilUpdater::update(double t)
  {
// get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();
// get input/output fields
    const Lucee::Field<2, double>& phi = this->getInp<Lucee::Field<2, double> >(0);
    const Lucee::Field<2, double>& aPar = this->getInp<Lucee::Field<2, double> >(1);
    Lucee::Field<2, double>& hamil = this->getOut<Lucee::Field<2, double> >(0);

    std::vector<int> ndIds;
    nodalBasis->getExclusiveNodeIndices(ndIds);
    unsigned numNodes = ndIds.size(); // number of nodes owned by each cell

    int idx[2];
    Lucee::Matrix<double> nodeCoords(nodalBasis->getNumNodes(), 3);

    Lucee::ConstFieldPtr<double> phiPtr = phi.createConstPtr();
    Lucee::ConstFieldPtr<double> aParPtr = aPar.createConstPtr();
    Lucee::FieldPtr<double> hamilPtr = hamil.createPtr();

    Lucee::Region<2, int> localExtRgn = phi.getExtRegion();
    Lucee::RowMajorSequencer<2> seq(localExtRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      phi.setPtr(phiPtr, idx);
      aPar.setPtr(aParPtr, idx);
      hamil.setPtr(hamilPtr, idx);

      nodalBasis->setIndex(idx);
      nodalBasis->getNodalCoordinates(nodeCoords);

      for (unsigned n=0; n<numNodes; ++n)
      {
        unsigned nn = ndIds[n]; // local node number we are working on
// H = 1/2m*(w-q*A)**2 + q*phi
        double t1 = (nodeCoords(nn,1)-charge*aParPtr[n]);
        hamilPtr[n] = 0.5/mass*t1*t1 + charge*phiPtr[n];
      }
    }

    return Lucee::UpdaterStatus();
  }

  void
  NonLinEmGke1dHamilUpdater::declareTypes()
  {
// takes two inputs (phi, Apar)
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
// returns one output, nonlinear Hamiltonian
    this->appendOutVarType(typeid(Lucee::Field<2, double>));    
  }
}
