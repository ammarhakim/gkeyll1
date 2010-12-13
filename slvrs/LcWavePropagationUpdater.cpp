/**
 * @file	LcWavePropagationUpdater.cpp
 *
 * @brief	Wave propagation solver.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// lucee includes
#include <LcDirSequencer.h>
#include <LcField.h>
#include <LcStructuredGridBase.h>
#include <LcWavePropagationUpdater.h>

// std includes
#include <iostream>
#include <string>
#include <vector>

namespace Lucee
{
  static const unsigned NO_LIMITER = 0;
  static const unsigned MINMOD_LIMITER = 1;
  static const unsigned SUPERBEE_LIMITER = 2;
  static const unsigned VAN_LEER_LIMITER = 3;
  static const unsigned MONOTONIZED_CENTERED_LIMITER = 4;
  static const unsigned BEAM_WARMING_LIMITER = 5;

// set ids for module system
  template <> const char *WavePropagationUpdater<1>::id = "WavePropagation1D";
  template <> const char *WavePropagationUpdater<2>::id = "WavePropagation2D";
  template <> const char *WavePropagationUpdater<3>::id = "WavePropagation3D";

  template <unsigned NDIM>
  void
  WavePropagationUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    UpdaterIfc::readInput(tbl);
// equation to solve
    if (tbl.hasObject<Lucee::HyperEquation>("equation"))
      equation = &tbl.getObjectAsBase<Lucee::HyperEquation>("equation");
    else
    {
      Lucee::Except lce("WavePropagationUpdater::readInput: Must specify an equation to solve!");
      throw lce;
    }
    
// limiter to use
    limiter = NO_LIMITER; // by default no limiter
    if (tbl.hasString("limiter"))
    {
      std::string lim = tbl.getString("limiter");
      if (lim == "no-limiter")
        limiter = NO_LIMITER;
      else if (lim == "min-mod")
        limiter = MINMOD_LIMITER;
      else if (lim == "superbee")
        limiter = SUPERBEE_LIMITER;
      else if (lim == "van-leer")
        limiter = VAN_LEER_LIMITER;
      else if (lim == "monotonized-centered")
        limiter = MONOTONIZED_CENTERED_LIMITER;
      else if (lim == "beam-warming")
        limiter = BEAM_WARMING_LIMITER;
      else
      {
        Lucee::Except lce("WavePropagationUpdater::readInput: Do not recognize limiter type '");
        lce << lim << "'" << std::endl;
        throw lce;
      }
    }
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  WavePropagationUpdater<NDIM>::update(double t)
  {
// get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
// get input/output arrays
    const Lucee::Field<NDIM, double>& q = this->getInp<Lucee::Field<NDIM, double> >(0);
    Lucee::Field<NDIM, double>& qNew = this->getOut<Lucee::Field<NDIM, double> >(0);

// cell spacing
    double dx[NDIM];
    for (unsigned n=0; n<NDIM; ++n) dx[n] = grid.getDx(n);
// time-step
    double dt = t-this->getCurrTime();
// local region to index
    Lucee::Region<NDIM, int> localRgn = grid.getLocalBox();

    unsigned meqn = equation->getNumEqns();
    unsigned mwave = equation->getNumWaves();

// indices
    int idx[NDIM], idxl[NDIM];
// pointers to data
    Lucee::ConstFieldPtr<double> qPtr = q.createConstPtr();
    Lucee::ConstFieldPtr<double> qPtrl = q.createConstPtr();
    Lucee::FieldPtr<double> qNewPtr = qNew.createPtr();
    Lucee::FieldPtr<double> qNewPtrl = qNew.createPtr();
// jumps
    Lucee::FieldPtr<double> jump(meqn);
// speeds
    Lucee::FieldPtr<double> s(mwave);
// waves
    Lucee::Matrix<double> waves(meqn, mwave);

// loop, updating slices in each dimension
    for (unsigned dir=0; dir<NDIM; ++dir)
    {
// create sequencer for looping over slices in 'dir' direction. Also
// make stencil with one cell on each side
      Lucee::DirSequencer<NDIM> seq(localRgn, dir, 1, 1);
      
      while (seq.step())
      {
// get index of this cell
        seq.fillWithIndex(0, idx);
// get index of left cell
        seq.fillWithIndex(-1, idxl);

// get hold of solution in these cells
        q.setPtr(qPtr, idx);
        q.setPtr(qPtrl, idxl);
// compute jump
        for (unsigned m=0; m<meqn; ++m)
          jump[m] = qPtr[m] - qPtrl[m];

// calculate waves and speeds
        equation->waves(jump, qPtrl, qPtr, waves, s);
// compute fluctuations


      }
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  WavePropagationUpdater<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class WavePropagationUpdater<1>;
  template class WavePropagationUpdater<2>;
  template class WavePropagationUpdater<3>;
}
