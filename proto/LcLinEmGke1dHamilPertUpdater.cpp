/**
 * @file	LcLinEmGke1dPertHamilUpdater.cpp
 *
 * @brief	Compute linearized perturbed Hamiltonian
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLinEmGke1dHamilPertUpdater.h>

namespace Lucee
{
  const char *LinEmGke1dPertHamilUpdater::id = "LinEmGke1DPertHamil";

  LinEmGke1dPertHamilUpdater::LinEmGke1dPertHamilUpdater()
    : UpdaterIfc()
  {
  }

  void
  LinEmGke1dPertHamilUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);
    charge = tbl.getNumber("charge");
    mass = tbl.getNumber("mass");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis");
    else
      throw Lucee::Except("EvalOnNodesUpdater::readInput: Must specify element to use using 'basis'");
  }

  void
  LinEmGke1dPertHamilUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus
  LinEmGke1dPertHamilUpdater::update(double t)
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
// H = -q/m*w*Apar + q*phi
        hamilPtr[n] = -charge/mass*nodeCoords(nn,1)*aParPtr[n] + charge*phiPtr[n];
      }
    }

    return Lucee::UpdaterStatus();
  }

  void
  LinEmGke1dPertHamilUpdater::declareTypes()
  {
// takes two inputs (phi, Apar)
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
// returns one output, perturbed Hamiltonian
    this->appendOutVarType(typeid(Lucee::Field<2, double>));    
  }
}
