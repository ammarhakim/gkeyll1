/**
 * @file	LcRadialStairSteppedBcUpdater.cpp
 *
 * @brief	Base class for boundary conditions for radial stair-stepped boundaries.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGlobals.h>
#include <LcAlignedRectCoordSys.h>
#include <LcRadialStairSteppedBcUpdater.h>
#include <LcField.h>
#include <LcStructuredGridBase.h>

// txbase includes
#include <TxHdf5Base.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
// right/left flags
  enum {RSSB_REAL = -1, RSSB_SKINGHOST = 0, RSSB_MIDDLE = 1};

  template <> const char *RadialStairSteppedBcUpdater<1>::id = "RadialStairSteppedBc1D";
  template <> const char *RadialStairSteppedBcUpdater<2>::id = "RadialStairSteppedBc2D";
  template <> const char *RadialStairSteppedBcUpdater<3>::id = "RadialStairSteppedBc3D";

  template <unsigned NDIM>
  RadialStairSteppedBcUpdater<NDIM>::RadialStairSteppedBcUpdater() 
    : Lucee::UpdaterIfc() 
  {
  }

  template <unsigned NDIM>
  RadialStairSteppedBcUpdater<NDIM>::~RadialStairSteppedBcUpdater() 
  {
    delete rssBnd;
    isSafeToWrite = false;
  }

  template <unsigned NDIM>
  void
  RadialStairSteppedBcUpdater<NDIM>::readInput(Lucee::LuaTable& tbl) 
  {
    UpdaterIfc::readInput(tbl);

    bcDir = 0;
// read direction to apply
    if (tbl.hasNumber("dir"))
      bcDir = tbl.getNumber("dir");

    origin = 0.;
    if (tbl.hasNumVec("origin"))
    {
      std::vector<double> originVec = tbl.getNumVec("origin");
      if (originVec.size() != NDIM)
      {
        Lucee::Except lce("RadialStairSteppedBcUpdater::readInput: 'origin' should have exactly");
        lce << NDIM << " elements. Instead has " << originVec.size() << std::endl;
        throw lce;
      }
      for (unsigned d = 0; d < NDIM; ++d)
      {
        origin[d] = originVec[d];
      }
    }

    if (!tbl.hasNumber("radius"))
      Lucee::Except lce("RadialStairSteppedBcUpdater::readInput: Must specify 'radius' of stair stepped boundary!");

    radius = tbl.getNumber("radius");
    if (radius <= 0)
      Lucee::Except lce("RadialStairSteppedBcUpdater::readInput: 'radius' of stair stepped boundary must be greater than zero!");

    if (tbl.hasNumVec("constComponents"))
    {
// get list of components to apply this boundary condition
      std::vector<double> cd = tbl.getNumVec("constyComponents");
      for (unsigned i=0; i<cd.size(); ++i)
        constComponents.push_back( (int) cd[i] );
    }

    if (constComponents.size() > 0)
    {
// read constant values to set to
      if (tbl.hasNumVec("constValues"))
      {
        constValues = tbl.getNumVec("constValues");
        if (constValues.size() != constComponents.size())
          throw Lucee::Except(
            "RadialStairSteppedBcUpdater::readInput: 'constValues' table must have same size as 'constComponents' table");
      }
      else
      {
        throw Lucee::Except("RadialStairSteppedBcUpdater::readInput: 'constValues' table must be specified when constComponents is not empty");
      }
    }

    if (tbl.hasNumVec("copyComponents"))
    {
// get list of components to apply this boundary condition
      std::vector<double> cd = tbl.getNumVec("copyComponents");
      for (unsigned i=0; i<cd.size(); ++i)
        copyComponents.push_back( (int) cd[i] );
    }

    if (tbl.hasNumVec("reflectComponents"))
    {
// get list of first components of vectors to apply this boundary condition
      std::vector<double> cd = tbl.getNumVec("reflectComponents");
      for (unsigned i=0; i<cd.size(); ++i)
        reflectComponents.push_back( (int) cd[i] );
    }

    if (tbl.hasNumVec("absorbComponents"))
    {
// get list of first components of vectors to apply this boundary condition
      std::vector<double> cd = tbl.getNumVec("absorbComponents");
      for (unsigned i=0; i<cd.size(); ++i)
        absorbComponents.push_back( (int) cd[i] );
    }

    if (tbl.hasNumVec("reflectTransComponents"))
    {
// get list of first components of vectors to apply this boundary condition
      std::vector<double> cd = tbl.getNumVec("reflectTransComponents");
      for (unsigned i=0; i<cd.size(); ++i)
        reflectTransComponents.push_back( (int) cd[i] );
    }

    if (tbl.hasNumVec("zeroRadialComponents"))
    {
// get list of first components of vectors to apply this boundary condition
      std::vector<double> cd = tbl.getNumVec("zeroRadialComponents");
      for (unsigned i=0; i<cd.size(); ++i)
        zeroRadialComponents.push_back( (int) cd[i] );
    }

    if (tbl.hasNumVec("zeroTransComponents"))
    {
// get list of first components of vectors to apply this boundary condition
      std::vector<double> cd = tbl.getNumVec("zeroTransComponents");
      for (unsigned i=0; i<cd.size(); ++i)
        zeroTransComponents.push_back( (int) cd[i] );
    }

// set pointer to in/out field
    inOut = &tbl.getObject<Lucee::Field<NDIM, double> >("inOutField");
  }

  template <unsigned NDIM>
  void
  RadialStairSteppedBcUpdater<NDIM>::initialize() 
  {
// call base class method
    UpdaterIfc::initialize();

    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    int lg[NDIM], ug[NDIM];
    for (unsigned d=0; d<NDIM; ++d)
      lg[d] = ug[d] = 2;
    rssBnd = new Lucee::Field<NDIM, double>(localRgn, 1, lg, ug); // to store boundary information
    isSafeToWrite = false;

    Lucee::FieldPtr<double> iop = inOut->createPtr();
    Lucee::FieldPtr<double> iopm = inOut->createPtr();

    (*rssBnd) = RSSB_MIDDLE; // clear out field

    Lucee::FieldPtr<double> rssp = rssBnd->createPtr();

// loop, figuring out where and how the boundary edge is oriented
    for (unsigned dir=0; dir<NDIM; ++dir)
    {
// create sequencer to loop over *each* 1D slice in 'dir' direction
      Lucee::RowMajorSequencer<NDIM> seq(localRgn.deflate(dir));

// lower and upper bounds of 1D slice
      int sliceLower = localRgn.getLower(dir)-1;
      int sliceUpper = localRgn.getUpper(dir)+1;

      while (seq.step())
      {
        int idx[NDIM], idxm[NDIM];
        seq.fillWithIndex(idx);
        seq.fillWithIndex(idxm);
// loop over each slice
        for (int i=sliceLower; i<sliceUpper; ++i)
        {
          idx[dir] = i;
          inOut->setPtr(iop, idx);
          rssBnd->setPtr(rssp, idx);

          if (iop[0] >0)
          {
// if current cell is real, mark it and skip to next cell
            rssp[0] = RSSB_REAL;
            continue;
          }

// now find skin ghost cells by finding if any adjacent cell is inside domain
// lower and upper bounds of adjacent cells
          int adjLower[NDIM], adjUpper[NDIM];
          for (unsigned d = 0; d < NDIM; ++d)
          {
            adjLower[d] = idx[d] - 1;
            adjUpper[d] = idx[d] + 2;
          }
          Lucee::Region<NDIM, int> adjRgn(adjLower, adjUpper);
          Lucee::RowMajorSequencer<NDIM> adjSeq(adjRgn);

          bool foundSkin = false;
          while (adjSeq.step() && !foundSkin)
          {
            adjSeq.fillWithIndex(idxm);
            inOut->setPtr(iopm, idxm);
// if any one of the adjacent cells is inside the domain, the current cell is a skin ghost cell
            if (iopm[0] > 0)
            {
              rssp[0] = RSSB_SKINGHOST;
              foundSkin = true;
              break;
            }
          } // while adjSeq.step

        } // for i sliceLower to sliceUpper
      } // while seq.step
    } // for dir dim
    isSafeToWrite = true;
  }

  template <unsigned NDIM>
  int
  RadialStairSteppedBcUpdater<NDIM>::luaWrite(lua_State *L)
  {
    RadialStairSteppedBcUpdater<NDIM> *updater
      = Lucee::PointerHolder<RadialStairSteppedBcUpdater<NDIM> >::getObj(L);
    std::string nm = lua_tostring(L, 2);
    if (updater->isSafeToWrite)
      updater->write(nm);

    return 0;
  }

  template <unsigned NDIM>
  void
  RadialStairSteppedBcUpdater<NDIM>::write(const std::string& nm)
  { 
    std::string outPrefix = Loki::SingletonHolder<Lucee::Globals>::Instance().outPrefix;
    std::string outNm = outPrefix + "_" + nm;
    TxCommBase& comm = *this->getComm();

    TxIoBase *io = new TxHdf5Base(&comm);
    TxIoNodeType fn = io->createFile(outNm);
    rssBnd->writeToFile(*io, fn, this->getName());
    delete io;
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  RadialStairSteppedBcUpdater<NDIM>::update(double t) 
  {
// don't do anything if apply direction is not in range
    if (bcDir>=NDIM)
    {
      Lucee::Except lce("RadialStairSteppedBcUpdater::update: Update direction should be less than ");
      lce << NDIM << ". Provided " << bcDir << " instead";
      throw lce;
    }

    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
// get field to update
    Lucee::Field<NDIM, double>& A = this->getOut<Lucee::Field<NDIM, double> >(0);
    unsigned numComponents = A.getNumComponents();

// create coordinate system along this direction (note that this
// updater only works on Cartesian meshes. Hence, using an aligned CS
// is perfectly fine)
    Lucee::AlignedRectCoordSys coordSys(bcDir);

    Lucee::FieldPtr<double> Aptr = A.createPtr();
    Lucee::FieldPtr<double> Aptrm = A.createPtr();
    Lucee::ConstFieldPtr<double> rssp = rssBnd->createConstPtr();
    Lucee::ConstFieldPtr<double> rsspm = rssBnd->createConstPtr();
    double xc[3], xcm[3];

// create sequencer to loop over *each* 1D slice in 'dir' direction
    Lucee::RowMajorSequencer<NDIM> seq(localRgn.deflate(bcDir));

// lower and upper bounds of 1D slice
    int sliceLower = localRgn.getLower(bcDir)-1;
    int sliceUpper = localRgn.getUpper(bcDir)+1;

    while (seq.step())
    {
      int idx[NDIM], idxm[NDIM];
      seq.fillWithIndex(idx);
// loop over each slice
      for (int i=sliceLower; i<sliceUpper; ++i)
      {
// set current cell index
        idx[bcDir] = i;
        rssBnd->setPtr(rssp, idx);
          
        if (rssp[0] == RSSB_SKINGHOST)
        { // found a skin ghost cell
// get current centroid coordinates
          grid.setIndex(idx);
          grid.getCentroid(xc);
          A.setPtr(Aptr, idx);

// compute coordinate relative to origin
          for (unsigned d=0; d<NDIM; ++d)
            xc[d] -= origin[d];
// compute radius square of the current cell
          double rc2 = 0.;
          for (unsigned d=0; d<NDIM; ++d)
            rc2 += xc[d]*xc[d];

// lower and upper bounds of adjacent cells
          int adjLower[NDIM], adjUpper[NDIM];
          for (unsigned d = 0; d < NDIM; ++d)
          {
            adjLower[d] = idx[d] - 1;
            adjUpper[d] = idx[d] + 2;
          }
          Lucee::Region<NDIM, int> adjRgn(adjLower, adjUpper);
          Lucee::RowMajorSequencer<NDIM> adjSeq(adjRgn);

          unsigned cnt = 0; // count of adjacent cells
          double wt; // weight of each cell
          double sum = 0.; // sum of weights
          std::vector<double> qFromAll (numComponents, 0.); // container for contribution from all cells
          while (adjSeq.step())
          {
// set adjacent cell index
            adjSeq.fillWithIndex(idxm);
// set pointers for adjacent cell
            rssBnd->setPtr(rsspm, idxm);

            if (rsspm[0] == RSSB_REAL)
            { // found an adjacent cell inside the domain
// get centroid coordinate for adjacent cell
              grid.setIndex(idxm);
              grid.getCentroid(xcm);
// set pointers for adjacent cell
              A.setPtr(Aptrm, idxm);

// compute coordinate relative to origin
              for (unsigned d=0; d<NDIM; ++d)
                xcm[d] -= origin[d];

// use equal weights for now; might depend on xcm
              wt = 1.;
              sum = sum + wt;
              for (unsigned i = 0; i < numComponents; ++i)
// accumulate contributions from all adjacent cells inside the domain
                qFromAll[i] += Aptrm[i] * wt;
            }
            cnt = cnt + 1;
          }
          for (unsigned i = 0; i < numComponents; ++i)
// compute the weighted average
            qFromAll[i] /= sum;

          for (unsigned i=0; i<this->constComponents.size(); ++i)
          {
// simply const
              Aptr[constComponents[i]] = constValues[i];
          }

          for (unsigned i=0; i<this->copyComponents.size(); ++i)
          {
// simply copy
            unsigned c = copyComponents[i];
            Aptr[c] = qFromAll[c];
          }
          
          for (unsigned i=0; i<this->reflectComponents.size(); ++i)
          {
// copy the transverse component and reflect the radial component
            unsigned c = reflectComponents[i];
            double v_radial_r = 0.; // v_radial*rc
            for (unsigned d=0; d<NDIM; ++d)
            {
              v_radial_r += qFromAll[c+d]*xc[d];
            }
          
            for (unsigned d=0; d<NDIM; ++d)
            {
// compute cartesian components of the reflected vector according to
// v[i] = v[i] - 2*v_r*x[i]/r
              Aptr[c+d] = qFromAll[c+d] - 2.*v_radial_r*xc[d]/rc2;
            }
          }
          
          for (unsigned i=0; i<this->absorbComponents.size(); ++i)
          {
// copy the transverse component and reflect the radial component when it is inward
            unsigned c = absorbComponents[i];
            double v_radial_r = 0.; // v_radial*rc
            for (unsigned d=0; d<NDIM; ++d)
            {
              v_radial_r += qFromAll[c+d]*xc[d];
            }
            if (v_radial_r > 0)
// reflect the vector if radial component is outward
              for (unsigned d=0; d<NDIM; ++d)
              {
// compute cccartesian components of the reflected vector according to
// v[i] = v[[[i] - 2*v_r*x[i]/r
                Aptr[c+d] = qFromAll[c+d] - 2.*v_radial_r*xc[d]/rc2;
              }
            else
// copy the vector if radial component is inward
              for (unsigned d=0; d<NDIM; ++d)
              {
                Aptr[c+d] = qFromAll[c+d];
              }
          }
             
          for (unsigned i=0; i<this->reflectTransComponents.size(); ++i)
          {
// copy the radial component and reflect the transverse component
            unsigned c = reflectTransComponents[i];
            double v_radial_r = 0.; // v_radial*rc
            for (unsigned d=0; d<NDIM; ++d)
            {
              v_radial_r += qFromAll[c+d]*xc[d];
            }
          
            for (unsigned d=0; d<NDIM; ++d)
            {
// compute cartesian components of the reflected vector according to
// v[i] = - v[i] + 2*v_r*x[i]/r
              Aptr[c+d] = - qFromAll[c+d] + 2.*v_radial_r*xc[d]/rc2;
            }
          }

          for (unsigned i=0; i<this->zeroRadialComponents.size(); ++i)
          {
// copy the transverse component zero radial component
            unsigned c = zeroRadialComponents[i];
            double v_radial_r = 0.; // v_radial*rc
            for (unsigned d=0; d<NDIM; ++d)
            {
              v_radial_r += qFromAll[c+d]*xc[d];
            }
          
            for (unsigned d=0; d<NDIM; ++d)
            {
// compute cartesian components of the reflected vector according to
// v[i] = v[i] - v_r*x[i]/r
              Aptr[c+d] = qFromAll[c+d] - v_radial_r*xc[d]/rc2;
            }
          }
        
          for (unsigned i=0; i<this->zeroTransComponents.size(); ++i)
          {
// copy the radial component with zero transverse component
            unsigned c = zeroTransComponents[i];
            double v_radial_r = 0.; // v_radial*rc
            for (unsigned d=0; d<NDIM; ++d)
            {
              v_radial_r += qFromAll[c+d]*xc[d];
            }
          
            for (unsigned d=0; d<NDIM; ++d)
            {
// compute cartesian components of the reflected vector according to
// v[i] = v_r*x[i]/r
              Aptr[c+d] = v_radial_r*xc[d]/rc2;
            }
          }

        } // if current cell is skin ghost
      }
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  RadialStairSteppedBcUpdater<NDIM>::declareTypes() 
  {
// any number of output fields
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void
  RadialStairSteppedBcUpdater<NDIM>::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
// call base class Lua methods
    UpdaterIfc::appendLuaCallableMethods(lfm);

    lfm.appendFunc("setDir", luaSetDir);
    lfm.appendFunc("write", luaWrite);
  }

  template <unsigned NDIM>
  int
  RadialStairSteppedBcUpdater<NDIM>::luaSetDir(lua_State *L)
  {
    RadialStairSteppedBcUpdater<NDIM> *updater
      = Lucee::PointerHolder<RadialStairSteppedBcUpdater<NDIM> >::getObj(L);
    int d = (unsigned) lua_tonumber(L, 2); // current time to set
    updater->setDir(d);

    return 0;
  }

// instantiations
  template class RadialStairSteppedBcUpdater<1>;
  template class RadialStairSteppedBcUpdater<2>;
  template class RadialStairSteppedBcUpdater<3>;
}
