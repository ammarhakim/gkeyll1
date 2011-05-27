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
#include <LcAlignedRectCoordSys.h>
#include <LcDirSequencer.h>
#include <LcField.h>
#include <LcMathLib.h>
#include <LcStructuredGridBase.h>
#include <LcWavePropagationUpdater.h>

// std includes
#include <algorithm>
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
  static const unsigned ZERO_LIMITER = 6;

// set ids for module system
  template <> const char *WavePropagationUpdater<1>::id = "WavePropagation1D";
  template <> const char *WavePropagationUpdater<2>::id = "WavePropagation2D";
  template <> const char *WavePropagationUpdater<3>::id = "WavePropagation3D";

  template <unsigned NDIM>
  void
  WavePropagationUpdater<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

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

// directions to update
    if (tbl.hasNumVec("updateDirections"))
    {
      std::vector<double> ud = tbl.getNumVec("updateDirections");
      for (unsigned i=0; i<std::min<unsigned>(3, ud.size()); ++i)
      {
        unsigned d = (unsigned) ud[i];
        if (d<3)
          updateDims.push_back(d);
        else
          throw Lucee::Except("updateDirections must be a table with 0, 1, or 2");
      }
    }
    else
    {
      for (unsigned i=0; i<NDIM; ++i)
        updateDims.push_back(i);
    }

    cfl = tbl.getNumber("cfl"); // CFL number
    cflm = tbl.getNumber("cflm"); // maximum CFL number
    
// limiter to use
    limiter = NO_LIMITER; // by default no limiter
    if (tbl.hasString("limiter"))
    {
      std::string lim = tbl.getString("limiter");
      if (lim == "no-limiter")
        limiter = NO_LIMITER;
      else if (lim == "minmod")
        limiter = MINMOD_LIMITER;
      else if (lim == "superbee")
        limiter = SUPERBEE_LIMITER;
      else if (lim == "van-leer")
        limiter = VAN_LEER_LIMITER;
      else if (lim == "monotonized-centered")
        limiter = MONOTONIZED_CENTERED_LIMITER;
      else if (lim == "beam-warming")
        limiter = BEAM_WARMING_LIMITER;
      else if (lim == "zero")
        limiter = ZERO_LIMITER;
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

// ghost cells along slice
    int lg[1], ug[1];
    lg[0] = 2; ug[0] = 2;

// allocate data
    for (unsigned dir=0; dir<NDIM; ++dir)
    {
      int lower[1], upper[1];
      lower[0] = localRgn.getLower(dir); 
      upper[0] = localRgn.getUpper(dir);
      Lucee::Region<1, int> slice(lower, upper);
      apdq.push_back(Lucee::Field<1, double>(slice, meqn, lg, ug));
      amdq.push_back(Lucee::Field<1, double>(slice, meqn, lg, ug));
      speeds.push_back(Lucee::Field<1, double>(slice, mwave, lg, ug));
      waves.push_back(Lucee::Field<1, double>(slice, meqn*mwave, lg, ug));
      fs.push_back(Lucee::Field<1, double>(slice, meqn, lg, ug));
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

// time-step
    double dt = t-this->getCurrTime();
// local region to index
    Lucee::Region<NDIM, int> localRgn = grid.getLocalBox();

    unsigned meqn = equation->getNumEqns();
    unsigned mwave = equation->getNumWaves();

// pointers to data
    Lucee::ConstFieldPtr<double> qPtr = q.createConstPtr();
    Lucee::ConstFieldPtr<double> qPtrl = q.createConstPtr();
    Lucee::FieldPtr<double> qNewPtr = qNew.createPtr();
    Lucee::FieldPtr<double> qNewPtrl = qNew.createPtr();
// rotated left/right states
    Lucee::FieldPtr<double> qLocal(meqn), qLocall(meqn);
// waves in local coordinate system
    Lucee::Matrix<double> wavesLocal(meqn, mwave);

// to store jump across interface
    Lucee::FieldPtr<double> jump(meqn);

// maximum CFL number used
    double cfla = 0.0;

// loop, updating slices in each requested dimension
    for (unsigned d=0; d<updateDims.size(); ++d)
    {
      unsigned dir = updateDims[d]; // direction to update
      double dtdx = dt/grid.getDx(dir);
// create coordinate system along this direction
      Lucee::AlignedRectCoordSys coordSys(dir);

// create sequencer to loop over *each* 1D slice in 'dir' direction
      Lucee::RowMajorSequencer<NDIM> seq(localRgn.deflate(dir));

// pointers to data
      Lucee::FieldPtr<double> apdqPtr = apdq[dir].createPtr();
      Lucee::FieldPtr<double> amdqPtr = amdq[dir].createPtr();
      Lucee::FieldPtr<double> speedsPtr = speeds[dir].createPtr();
      Lucee::FieldPtr<double> wavesPtr = waves[dir].createPtr();
      Lucee::FieldPtr<double> fsPtr = fs[dir].createPtr();
      Lucee::FieldPtr<double> fsPtr1 = fs[dir].createPtr();

// lower and upper bounds of 1D slice. (We need to make sure that the
// Riemann problem is computed for one edge outside the domain
// interior. This is needed to limit the waves on the domain boundary)
      int sliceLower = localRgn.getLower(dir)-1;
      int sliceUpper = localRgn.getUpper(dir)+2;

// loop over each 1D slice
      while (seq.step())
      {
        int idx[NDIM], idxl[NDIM];
        seq.fillWithIndex(idx);
        seq.fillWithIndex(idxl);
// loop over each edge in slice
        for (int i=sliceLower; i<sliceUpper; ++i)
        {
          idx[dir] = i; // cell right of edge
          idxl[dir] = i-1; // cell left of edge
// get hold of solution in these cells
          q.setPtr(qPtr, idx);
          q.setPtr(qPtrl, idxl);
// rotate data to local coordinate on cell face
          equation->rotateToLocal(coordSys, &qPtr[0], &qLocal[0]);
          equation->rotateToLocal(coordSys, &qPtrl[0], &qLocall[0]);

// attach pointers to fluctuations, speeds, waves (note these are 1D arrays)
          apdq[dir].setPtr(apdqPtr, i);
          amdq[dir].setPtr(amdqPtr, i);
          speeds[dir].setPtr(speedsPtr, i);
          waves[dir].setPtr(wavesPtr, i);

// compute jump across edge
          for (unsigned m=0; m<meqn; ++m)
            jump[m] = qLocal[m] - qLocall[m];

// calculate waves and speeds
          equation->waves(coordSys, jump, qLocal, qLocall, wavesLocal, speedsPtr);
// rotate waves back to global frame (stored in waves[dir] array)
          Lucee::Matrix<double> wavesGlobal(meqn, mwave, wavesPtr);
          for (unsigned mw=0; mw<mwave; ++mw)
// the calls &wavesGlobal(0,mw) etc works because the columns (waves)
// are stored contiguously.
            equation->rotateToGlobal(coordSys, &wavesLocal(0,mw), &wavesGlobal(0,mw));

// compute fluctuations
          equation->qFluctuations(wavesGlobal, speedsPtr, amdqPtr, apdqPtr);

// compute first-order Gudonov update
          qNew.setPtr(qNewPtr, idx);
          qNew.setPtr(qNewPtrl, idxl);
          for (unsigned m=0; m<meqn; ++m)
          {
            qNewPtr[m] += -dtdx*apdqPtr[m];
            qNewPtrl[m] += -dtdx*amdqPtr[m];
          }

// compute CFL number used in this step
          for (unsigned mw=0; mw<mwave; ++mw)
            cfla = Lucee::max3(cfla, dtdx*speedsPtr[mw], -dtdx*speedsPtr[mw]);

        }
// check if time-step was too large
        if (cfla > cflm)
          return Lucee::UpdaterStatus(false, dt*cfl/cfla);

// apply limiters
        applyLimiters(waves[dir], speeds[dir]);

// compute second order corrections to flux (we need to go one cell
// beyond the last cell to ensure the right most edge flux is
// computed)
        for (int i=localRgn.getLower(dir); i<localRgn.getUpper(dir)+1; ++i)
        {
          fs[dir].setPtr(fsPtr, i);
          speeds[dir].setPtr(speedsPtr, i);
          waves[dir].setPtr(wavesPtr, i);

// create matrix to store waves
          Lucee::Matrix<double> wavesMat(meqn, mwave, wavesPtr);

          for (unsigned m=0; m<meqn; ++m)
          { // compute correction
            fsPtr[m] = 0.0;
            for (unsigned mw=0; mw<mwave; ++mw)
              fsPtr[m] += 0.5*std::abs(speedsPtr[mw])*(1.0 -
                std::abs(speedsPtr[mw])*dtdx)*wavesMat(m, mw);
          }
        }

// accumulate second order corrections
        for (int i=localRgn.getLower(dir); i<localRgn.getUpper(dir); ++i)
        {
          idx[dir] = i; // left edge of cell
// set pointers to proper locations
          qNew.setPtr(qNewPtr, idx);
          fs[dir].setPtr(fsPtr, i);
          fs[dir].setPtr(fsPtr1, i+1);

          for (unsigned m=0; m<meqn; ++m)
            qNewPtr[m] += -dtdx*(fsPtr1[m] - fsPtr[m]);
        }
      }
    }
    return Lucee::UpdaterStatus(true, dt*cfl/cfla);
  }

  template <unsigned NDIM>
  void
  WavePropagationUpdater<NDIM>::applyLimiters(
    Lucee::Field<1, double>& ws, const Lucee::Field<1, double>& sp)
  {
    double c, r, dotr, dotl, wnorm2, wlimitr=1;

    unsigned meqn = equation->getNumEqns();
    unsigned mwave = equation->getNumWaves();
// create indexer to access waves, stored as a meqn X mwave matrix
    int start[2] = {0, 0};
    unsigned shape[2] = {meqn, mwave};
    Lucee::ColMajorIndexer<2> idx(shape, start);

    Lucee::ConstFieldPtr<double> spPtr = sp.createConstPtr();

    int sliceLower = sp.getLower(0);
    int sliceUpper = sp.getUpper(0) + 1;
    for (unsigned mw=0; mw<mwave; ++mw)
    {
      dotr = 0.0;
// compute initial dotr value (this will become dotl in the loop)
      for (unsigned m=0; m<meqn; ++m)
        dotr += ws(sliceLower-1, idx.getIndex(m, mw))*ws(sliceLower, idx.getIndex(m, mw));

      for (int i=sliceLower; i<sliceUpper; ++i)
      {
        sp.setPtr(spPtr, i);
        wnorm2 = 0.0;
        dotl = dotr;
        dotr = 0.0;
        for (unsigned m=0; m<meqn; ++m)
        { // compute norm and dotr
          wnorm2 += ws(i, idx.getIndex(m, mw))*ws(i, idx.getIndex(m, mw));
          dotr += ws(i, idx.getIndex(m, mw))*ws(i+1, idx.getIndex(m, mw));
        }
        if (wnorm2 > 0.0)
        {
          if (spPtr[mw] > 0)
            r = dotl/wnorm2;
          else
            r = dotr/wnorm2;

          switch (limiter)
          {
            case NO_LIMITER:
                wlimitr = 1.0;
                break;
                
            case MINMOD_LIMITER:
                wlimitr = std::max(0.0, std::min(1.0, r));
                break;

            case SUPERBEE_LIMITER:
                wlimitr = Lucee::max3(0.0, std::min(1.0, 2*r), std::min(2.0, r));
                break;

            case VAN_LEER_LIMITER:
                wlimitr = (r+std::abs(r))/(1+std::abs(r));
                break;

            case MONOTONIZED_CENTERED_LIMITER:
                c = (1.+r)/2.;
                wlimitr = std::max(0.0, Lucee::min3(c, 2., 2.*r));
                break;

            case BEAM_WARMING_LIMITER:
                wlimitr = r;
                break;

            case ZERO_LIMITER:
                wlimitr = 0;
                break;

            default:
                ;
          }
// apply limiter
          for (unsigned m=0; m<meqn; ++m)
            ws(i, idx.getIndex(m, mw)) = wlimitr*ws(i, idx.getIndex(m, mw));
        }
      }
    }
  }

// instantiations
  template class WavePropagationUpdater<1>;
  template class WavePropagationUpdater<2>;
  template class WavePropagationUpdater<3>;
}
