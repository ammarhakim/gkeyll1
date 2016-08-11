/**
 * @file	LcSheathParticleSource1x1v.cpp
 *
 * @brief	Compute particle source to compensate for loss at boundaries
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGlobals.h>
#include <LcMathPhysConstants.h>
#include <LcSheathParticleSource1x1v.h>

// loki includes
#include <loki/Singleton.h>
#include <cmath>

// std includes
#include <vector>

namespace Lucee
{
  const char *SheathParticleSource1x1v::id = "SheathParticleSource1x1v";

  SheathParticleSource1x1v::SheathParticleSource1x1v()
  {
  }

  void
  SheathParticleSource1x1v::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis");
    else
      throw Lucee::Except("SheathParticleSource1x1v::readInput: Must specify element to use using 'basis'");

    vThermal = tbl.getNumber("thermalVelocity");
    uDrift = tbl.getNumber("driftVelocity");
  }

  void
  SheathParticleSource1x1v::initialize()
  {
    UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus
  SheathParticleSource1x1v::update(double t)
  {
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();
// get field with momentum fluxes
    const Lucee::Field<1, double>& ptclMom = this->getInp<Lucee::Field<1, double> >(0);
// get output source array
    Lucee::Field<2, double>& src = this->getOut<Lucee::Field<2, double> >(0);
    src = 0.0;

    unsigned numNodes = nodalBasis->getNumNodes();
    int idx[2];
    Lucee::Matrix<double> nodeCoords(nodalBasis->getNumNodes(), 3);
    Lucee::FieldPtr<double> ptr = src.createPtr();

    Lucee::Region<1, int> globalRgn = ptclMom.getGlobalRegion();
    Lucee::Region<2, double> compSpace = grid.getComputationalSpace();
    double Lx = compSpace.getShape(0);

    Lucee::ConstFieldPtr<double> ptclMomPtr = ptclMom.createConstPtr();
// compute net outflow of momentum
    ptclMom.setPtr(ptclMomPtr, globalRgn.getLower(0));
    double leftMomFlux = ptclMomPtr[0]; // left most node
    ptclMom.setPtr(ptclMomPtr, globalRgn.getUpper(0)-1);
    double rightMomFlux = ptclMomPtr[ptclMomPtr.getNumComponents()-1]; // right most node

// source balances outflow of momentum
    double Sn = (rightMomFlux-leftMomFlux)/Lx;
    double vt2 = vThermal*vThermal;
    double maxFact = Sn/std::sqrt(2*Lucee::PI*vt2);

    Lucee::Region<2, int> localExtRgn = src.getExtRegion();
    Lucee::RowMajorSequencer<2> seq(localExtRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      src.setPtr(ptr, idx);

      nodalBasis->setIndex(idx);
      nodalBasis->getNodalCoordinates(nodeCoords);

      for (unsigned n=0; n<numNodes; ++n)
      {
        double vu = nodeCoords(n,1)-uDrift;
        ptr[n] = maxFact*std::exp(-vu*vu/2/vt2);
      }
    }
    return Lucee::UpdaterStatus();
  }

  void
  SheathParticleSource1x1v::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }
}
