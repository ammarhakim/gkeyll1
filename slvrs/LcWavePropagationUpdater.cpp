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

    cfl = tbl.getNumber("cfl"); // CFL number
    cflm = tbl.getNumber("cfl"); // maximum CFL number
    
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
  void
  WavePropagationUpdater<NDIM>::initialize()
  {
// call base class method
    UpdaterIfc::initialize();
// get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
// local region to index
    Lucee::Region<NDIM, int> localRgn = grid.getLocalBox();
// determine number of equations and waves
    unsigned meqn = equation->getNumEqns();
    unsigned mwave = equation->getNumWaves();

// allocate data
    for (unsigned i=0; i<NDIM; ++i)
    {
      int lower[1], upper[1];
      lower[0] = localRgn.getLower(0); upper[0] = localRgn.getUpper(0);
      Lucee::Region<1, int> slice(lower, upper);
      apdq.push_back(Lucee::Field<1, double>(slice, meqn));
      amdq.push_back(Lucee::Field<1, double>(slice, meqn));
      speeds.push_back(Lucee::Field<1, double>(slice, mwave));
      waves.push_back(Lucee::Field<1, double>(slice, meqn*mwave));
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

// qnew <- qold
    qNew.copy(q);

// cell spacing
    double dx[NDIM];
    for (unsigned n=0; n<NDIM; ++n) 
      dx[n] = grid.getDx(n);
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

// to store jump across interface
    Lucee::FieldPtr<double> jump(meqn);

// loop, updating slices in each dimension
    for (unsigned dir=0; dir<NDIM; ++dir)
    {
// create sequencer to loop over *each* 1D slice in 'dir' direction
      Lucee::RowMajorSequencer<NDIM> seq(localRgn.deflate(dir));

// pointers to data
      Lucee::FieldPtr<double> apdqPtr = apdq[dir].createPtr();
      Lucee::FieldPtr<double> amdqPtr = amdq[dir].createPtr();
      Lucee::FieldPtr<double> speedsPtr = speeds[dir].createPtr();
      Lucee::FieldPtr<double> wavesPtr = waves[dir].createPtr();
      
// loop over each 1D slice
      while (seq.step())
      {
        int idx[NDIM], idxl[NDIM];
        seq.fillWithIndex(idx);
        seq.fillWithIndex(idxl);
// loop over slice
        for (int i=localRgn.getLower(dir); i<localRgn.getUpper(dir); ++i)
        {
          idx[dir] = i; // cell right of edge
          idxl[dir] = i-1; // cell left of edge
// get hold of solution in these cells
          q.setPtr(qPtr, idx);
          q.setPtr(qPtrl, idxl);
// attach pointers to fluctuations, speeds, waves
          apdq[dir].setPtr(apdqPtr, idx);
          amdq[dir].setPtr(amdqPtr, idx);
          speeds[dir].setPtr(speedsPtr, idx);
          waves[dir].setPtr(wavesPtr, idx);
// create matrix to store waves
          Lucee::Matrix<double> wavesMat(meqn, mwave, wavesPtr);

// compute jump across edge
          for (unsigned m=0; m<meqn; ++m)
            jump[m] = qPtr[m] - qPtrl[m];

// calculate waves and speeds
          equation->waves(jump, qPtrl, qPtr, wavesMat, speedsPtr);
// compute fluctuations
          equation->qFluctuations(wavesMat, speedsPtr, apdqPtr, amdqPtr);
        }
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
