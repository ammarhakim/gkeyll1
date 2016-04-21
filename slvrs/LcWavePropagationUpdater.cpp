/**
 * @file	LcWavePropagationUpdater.cpp
 *
 * @brief	Wave propagation solver.
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

  static bool isOutside(const Lucee::ConstFieldPtr<double>& p)
  { return p[0]<0; }
  
  static
  void clearFieldVect(std::vector<Lucee::Field<1, double>* >& fldVec)
  {
    std::vector<Lucee::Field<1, double>* >::iterator itr
      = fldVec.begin();
    for ( ; itr != fldVec.end(); ++itr)
      delete *itr;
    fldVec.clear();
  }

  template <unsigned NDIM>
  WavePropagationUpdater<NDIM>::~WavePropagationUpdater()
  {
    clearFieldVect(apdq);
    clearFieldVect(amdq);
    clearFieldVect(speeds);
    clearFieldVect(waves);
    clearFieldVect(fs);
  }

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
    cflm = 1.1*cfl; // use slightly large max CFL if not explicitly specified
    if (tbl.hasNumber("cflm"))
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
    hasLimiterField = false;
// check if there is an embedded boundary
    if (tbl.hasBool("hasLimiterField"))
      hasLimiterField = tbl.getBool("hasLimiterField");

// fetch pointer to location-based limiter
    if (hasLimiterField)
      limiterField = &tbl.getObject<Lucee::Field<NDIM, double> >("limiterField");

    bool hasFluxBc = false;
    for (unsigned d=0; d<NDIM; ++d)
    {
      hasLowerFluxBc[d] = hasUpperFluxBc[d] = false;
    }
// check if any direction has a flux BCs
    if (tbl.hasNumVec("lowerFluxDirs"))
    {
      hasFluxBc = true;
      std::vector<double> v = tbl.getNumVec("lowerFluxDirs");
      for (unsigned i=0; i<v.size(); ++i)
        hasLowerFluxBc[(int) v[i]] = true;
    }
    if (tbl.hasNumVec("upperFluxDirs"))
    {
      hasFluxBc = true;
      std::vector<double> v = tbl.getNumVec("upperFluxDirs");
      for (unsigned i=0; i<v.size(); ++i)
        hasUpperFluxBc[(int) v[i]] = true;
    }

    hasSsBnd = false;
// check if there is an embedded boundary
    if (tbl.hasBool("hasStairSteppedBoundary"))
      hasSsBnd = tbl.getBool("hasStairSteppedBoundary");

// fetch pointer to in/out field if there is an embedded boundary
    if (hasSsBnd)
      inOut = &tbl.getObject<Lucee::Field<NDIM, double> >("inOutField");

    zeroLimiterSsBnd = false;
    if (hasSsBnd && tbl.hasBool("zeroLimiterSsBnd"))
      zeroLimiterSsBnd = tbl.getBool("zeroLimiterSsBnd");

// fetch pointer to flux bc field (if any)
    if (hasFluxBc)
      fluxBc = &tbl.getObject<Lucee::Field<NDIM, double> >("boundaryFluxField");
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
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
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
      apdq.push_back(new Lucee::Field<1, double>(slice, meqn, lg, ug));
      amdq.push_back(new Lucee::Field<1, double>(slice, meqn, lg, ug));
      speeds.push_back(new Lucee::Field<1, double>(slice, mwave, lg, ug));
      waves.push_back(new Lucee::Field<1, double>(slice, meqn*mwave, lg, ug));
      fs.push_back(new Lucee::Field<1, double>(slice, meqn, lg, ug));
    }
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  WavePropagationUpdater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    const Lucee::Field<NDIM, double>& q = this->getInp<Lucee::Field<NDIM, double> >(0);
    Lucee::Field<NDIM, double>& qNew = this->getOut<Lucee::Field<NDIM, double> >(0);

    qNew.copy(q);

    double dt = t-this->getCurrTime();
// local region to index
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    
    unsigned meqn = equation->getNumEqns();
    unsigned mwave = equation->getNumWaves();

    Lucee::ConstFieldPtr<double> qPtr = q.createConstPtr();
    Lucee::ConstFieldPtr<double> qPtrl = q.createConstPtr();
    Lucee::FieldPtr<double> qNewPtr = qNew.createPtr();
    Lucee::FieldPtr<double> qNewPtrl = qNew.createPtr();
    Lucee::FieldPtr<double> qLocal(meqn), qLocall(meqn);
    Lucee::Matrix<double> wavesLocal(meqn, mwave);

// determine number of auxillary variables
    unsigned numAuxVars = this->getNumInpVars()-1; // first is always conserved variable
// store them
    std::vector<const Lucee::Field<NDIM, double>* > auxVars;
    for (unsigned i=0; i<numAuxVars; ++i)
      auxVars.push_back( &this->getInp<Lucee::Field<NDIM, double> >(i+1) );
// compute number of equations in each auxillary variable
    std::vector<unsigned> numAuxEqns;
    for (unsigned i=0; i<numAuxVars; ++i)
      numAuxEqns.push_back( auxVars[i]->getNumComponents() );

    std::vector<const double *> auxQl(numAuxVars), auxQr(numAuxVars);
    std::vector<const double *> inAuxQl(numAuxVars), inAuxQr(numAuxVars);    

    Lucee::FieldPtr<double> jump(meqn);

// these pointers are to the inOut and fluxBc fields, but as it can be
// NULL, I am using q to make their pointers. The reason this works is
// a "bug" in the code. I.e. pointers do not check if they are set to
// their parent fields. (AHH 7/13/2014)
    Lucee::ConstFieldPtr<double> ioPtr = q.createConstPtr();
    Lucee::ConstFieldPtr<double> ioPtr1 = q.createConstPtr();
    Lucee::ConstFieldPtr<double> fluxBcPtr = q.createConstPtr();

    double cfla = 0.0; // maximum CFL number used

// loop, updating slices in each requested dimension
    for (unsigned d=0; d<updateDims.size(); ++d)
    {
      unsigned dir = updateDims[d]; // direction to update
// create coordinate system along this direction
      Lucee::AlignedRectCoordSys coordSys(dir);

// create sequencer to loop over *each* 1D slice in 'dir' direction
      Lucee::RowMajorSequencer<NDIM> seq(localRgn.deflate(dir));

      Lucee::FieldPtr<double> apdqPtr = apdq[dir]->createPtr();
      Lucee::FieldPtr<double> amdqPtr = amdq[dir]->createPtr();
      Lucee::FieldPtr<double> speedsPtr = speeds[dir]->createPtr();
      Lucee::FieldPtr<double> wavesPtr = waves[dir]->createPtr();
      Lucee::FieldPtr<double> fsPtr = fs[dir]->createPtr();
      Lucee::FieldPtr<double> fsPtr1 = fs[dir]->createPtr();

// lower and upper bounds of 1D slice. (We need to make sure that the
// Riemann problem is computed for one edge outside the domain
// interior. This is needed to limit the waves on the domain boundary)
      int sliceLower = localRgn.getLower(dir)-1;
      int sliceUpper = localRgn.getUpper(dir)+2;

// adjust lower/upper bounds if flux boundaries are specified, and we
// are on the proper ranks
      if (hasLowerFluxBc[dir] && (q.getLower(dir) == q.getGlobalLower(dir)))
        sliceLower = localRgn.getLower(dir)+1;
      if (hasUpperFluxBc[dir] && (q.getUpper(dir) == q.getGlobalUpper(dir)))
        sliceUpper = localRgn.getUpper(dir);

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

          if (hasSsBnd)
          {
// if both cells attached to this edge are outside domain, skip it
            inOut->setPtr(ioPtr, idx);
            inOut->setPtr(ioPtr1, idxl);
            if (isOutside(ioPtr) && isOutside(ioPtr1))
              continue; // skip to next cell
          }

// get hold of solution in these cells
          q.setPtr(qPtr, idx);
          q.setPtr(qPtrl, idxl);
// rotate data to local coordinate on cell face
          equation->rotateToLocal(coordSys, &qPtr[0], &qLocal[0]);
          equation->rotateToLocal(coordSys, &qPtrl[0], &qLocall[0]);

// create input auxiliary variables list
          for (unsigned a=0; a<numAuxVars; ++a)
          {
            Lucee::ConstFieldPtr<double> aPtr = auxVars[a]->createConstPtr();
            auxVars[a]->setPtr(aPtr, idx);
            auxQr[a] = &aPtr[0];
            
            auxVars[a]->setPtr(aPtr, idxl);
            auxQl[a] = &aPtr[0];

// NOTE: we do not rotate auxiliary variables as the equation system
// should do the rotations if needed. This is perhaps inconsistent,
// but the HyperEquation class interface need not be cluttered with
// yet another set of rotation functions.
          }          

// attach pointers to fluctuations, speeds, waves (note these are 1D arrays)
          apdq[dir]->setPtr(apdqPtr, i);
          amdq[dir]->setPtr(amdqPtr, i);
          speeds[dir]->setPtr(speedsPtr, i);
          waves[dir]->setPtr(wavesPtr, i);

// compute jump across edge
          for (unsigned m=0; m<meqn; ++m)
            jump[m] = qLocal[m] - qLocall[m];

// calculate waves and speeds
          equation->waves(coordSys, jump, qLocall, qLocal, auxQl, auxQr, wavesLocal, speedsPtr);
// rotate waves back to global frame (stored in waves[dir] array)
          Lucee::Matrix<double> wavesGlobal(meqn, mwave, wavesPtr);
          for (unsigned mw=0; mw<mwave; ++mw)
// the calls &wavesGlobal(0,mw) etc works because the columns (waves)
// are stored contiguously.
            equation->rotateToGlobal(coordSys, &wavesLocal(0,mw), &wavesGlobal(0,mw));

// compute fluctuations
          equation->qFluctuations(coordSys, qLocall, qLocal,
            wavesGlobal, speedsPtr, amdqPtr, apdqPtr);

// get surface area of edge and volumes of two cells attached to edge
          grid.setIndex(idx); // cell right of edge
          double surfArea = grid.getSurfArea(dir);
          double areaVol = surfArea/grid.getVolume();

          grid.setIndex(idxl); // cell left of edge
          double areaVoll = surfArea/grid.getVolume();

// compute first-order Gudonov update
          qNew.setPtr(qNewPtr, idx);
          qNew.setPtr(qNewPtrl, idxl);
          for (unsigned m=0; m<meqn; ++m)
          {
            qNewPtr[m] += -dt*areaVol*apdqPtr[m]; // contribution to cell right of edge
            qNewPtrl[m] += -dt*areaVoll*amdqPtr[m]; // contribution to cell left of edge
          }

// compute CFL number used in this step
          for (unsigned mw=0; mw<mwave; ++mw)
            cfla = Lucee::max3(cfla, 
              std::fabs(dt*areaVol*speedsPtr[mw]), std::fabs(dt*areaVoll*speedsPtr[mw]));

        }
// check if time-step was too large
        if (cfla > cflm)
          return Lucee::UpdaterStatus(false, dt*cfl/cfla);

        applyLimiters(dir, idx, *waves[dir], *speeds[dir]);

// We need to go one cell beyond the last cell to ensure the right
// most edge flux is computed
      int sliceLower = localRgn.getLower(dir);
      int sliceUpper = localRgn.getUpper(dir)+1;
// adjust lower/upper bounds if flux boundaries are specified
      if (hasLowerFluxBc[dir] && (q.getLower(dir) == q.getGlobalLower(dir)))
        sliceLower = localRgn.getLower(dir)+1;
      if (hasUpperFluxBc[dir] && (q.getUpper(dir) == q.getGlobalUpper(dir)))
        sliceUpper = localRgn.getUpper(dir);

// compute second order corrections to flux This loop is over edges.
        for (int i=sliceLower; i<sliceUpper; ++i)
        {
          if (hasSsBnd)
          {
// if both cells attached to this edge are outside the domain, do not
// compute second order correction
            idx[dir] = i; // right cell
            inOut->setPtr(ioPtr, idx);
            idx[dir] = i-1; // left cell
            inOut->setPtr(ioPtr1, idx);
            if (isOutside(ioPtr) && isOutside(ioPtr1))
              continue; // skip to next cell
          }

          idx[dir] = i; // cell right of edge
          grid.setIndex(idx);
          double surfArea = grid.getSurfArea(dir);
          double cellVol = grid.getVolume();

          idx[dir] = i-1; // cell left of edge
          grid.setIndex(idx);
          cellVol += grid.getVolume();
          double areaVol = surfArea/(0.5*cellVol);

          fs[dir]->setPtr(fsPtr, i);
          speeds[dir]->setPtr(speedsPtr, i);
          waves[dir]->setPtr(wavesPtr, i);

          Lucee::Matrix<double> wavesMat(meqn, mwave, wavesPtr);

          for (unsigned m=0; m<meqn; ++m)
          { // compute correction
            fsPtr[m] = 0.0;
            for (unsigned mw=0; mw<mwave; ++mw)
              fsPtr[m] += 0.5*std::abs(speedsPtr[mw])*(1.0 -
                std::abs(speedsPtr[mw])*dt*areaVol)*wavesMat(m, mw);
          }
        }

// if we have flux BC, fill them into the "corrected" flux array, so
// that they are accounted for in following loop over cells
        if (hasLowerFluxBc[dir] && (q.getLower(dir) == q.getGlobalLower(dir)))
        {
          idx[dir] = q.getLower(dir);
          fluxBc->setPtr(fluxBcPtr, idx);
          fs[dir]->setPtr(fsPtr, q.getLower(dir));
          for (unsigned m=0; m<meqn; ++m)
            fsPtr[m] = fluxBcPtr[m];
        }
        if (hasUpperFluxBc[dir] && (q.getUpper(dir) == q.getGlobalUpper(dir)))
        {
          idx[dir] = q.getUpper(dir);
          fluxBc->setPtr(fluxBcPtr, idx);
          fs[dir]->setPtr(fsPtr, q.getUpper(dir));
          for (unsigned m=0; m<meqn; ++m)
            fsPtr[m] = fluxBcPtr[m];
        }

// accumulate second order corrections (this loop is over cells)
        for (int i=localRgn.getLower(dir); i<localRgn.getUpper(dir); ++i)
        {
          if (hasSsBnd)
          {
// if cell is outside domain, do not update solution
            idx[dir] = i; // cell index
            inOut->setPtr(ioPtr, idx);
            if (isOutside(ioPtr))
              continue; // skip to next cell
          }

          idx[dir] = i+1; //  right edge of cell
          grid.setIndex(idx);
          double surfArea1 = grid.getSurfArea(dir);

          idx[dir] = i; // left edge of cell
          grid.setIndex(idx);
          double surfArea = grid.getSurfArea(dir);
          double cellVol = grid.getVolume();

          qNew.setPtr(qNewPtr, idx);
          fs[dir]->setPtr(fsPtr, i);
          fs[dir]->setPtr(fsPtr1, i+1);

          for (unsigned m=0; m<meqn; ++m)
            qNewPtr[m] += -dt/cellVol*(surfArea1*fsPtr1[m] - surfArea*fsPtr[m]);
        }
      }
    }  
    return Lucee::UpdaterStatus(true, dt*cfl/cfla);
  }

  template <unsigned NDIM>
  void
  WavePropagationUpdater<NDIM>::applyLimiters(unsigned dir, int cellIdx[NDIM],
    Lucee::Field<1, double>& ws, const Lucee::Field<1, double>& sp)
  {
    if (limiter == NO_LIMITER and !limiterField) return;

    double c, r, dotr, dotl, wnorm2, wlimitr=1;

    unsigned meqn = equation->getNumEqns();
    unsigned mwave = equation->getNumWaves();
// create indexer to access waves, stored as a meqn X mwave matrix
    int start[2] = {0, 0};
    unsigned shape[2] = {meqn, mwave};
    Lucee::ColMajorIndexer<2> idx(shape, start);

    Lucee::ConstFieldPtr<double> spPtr = sp.createConstPtr();

// these pointers are to the inOut field, but as it can be NULL, I am
// using q to make their pointers. The reason this works is a "bug" in
// the code. I.e. pointers do not check if they are set to their
// parent fields. (AHH 7/13/2014)
    Lucee::ConstFieldPtr<double> ioPtr = sp.createConstPtr();
    Lucee::ConstFieldPtr<double> ioPtr1 = sp.createConstPtr();
    Lucee::ConstFieldPtr<double> lmtPtr = sp.createConstPtr();

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

        unsigned myLimiter = limiter;
        if (hasLimiterField && limiterField)
        {
          cellIdx[dir] = i;
          limiterField->setPtr(lmtPtr, cellIdx);
          myLimiter = lmtPtr[0];
        }

        if (hasSsBnd)
        {
// if both cells attached to this edge are outside the domain, do not
// limit wave
//
// (I no longer recall why the following lines have been commented
// out. AHH July 2015)
           cellIdx[dir] = i; // right cell
           inOut->setPtr(ioPtr, cellIdx);
           cellIdx[dir] = i-1; // left cell
           inOut->setPtr(ioPtr1, cellIdx);
           if (isOutside(ioPtr) && isOutside(ioPtr1))
           {
          //   continue; // skip to next cell
           }
           else if (zeroLimiterSsBnd && (isOutside(ioPtr) || isOutside(ioPtr1)))
             myLimiter = ZERO_LIMITER;
        }

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

          switch (myLimiter)
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
