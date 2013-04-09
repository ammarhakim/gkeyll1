/**
 * @file	LcPointSourceIfc.cpp
 *
 * @brief	Interface to "point" sources.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcPointSourceIfc.h>

namespace Lucee
{
// set module name
  const char *PointSourceIfc::id = "PointSource";

  PointSourceIfc::PointSourceIfc(unsigned nInp, unsigned nOut, bool allowArb)
    : Lucee::BasicObj("PointSource"), nInp(nInp), nOut(nOut), allowArb(allowArb),
      inpComponents(nInp), outComponents(nOut), srcJac(nOut, nOut)
  {
    for (unsigned i=0; i<nInp; ++i)
      inpComponents[i] = i;
    for (unsigned i=0; i<nOut; ++i)
      outComponents[i] = i;
  }

  PointSourceIfc::~PointSourceIfc()
  {
  }

  void
  PointSourceIfc::readInput(Lucee::LuaTable& tbl)
  {
    if (tbl.hasNumVec("inpComponents"))
    {
// get input component list
      std::vector<double> inpListDbl = tbl.getNumVec("inpComponents");
      if (allowArb == false)
      {
        if (inpListDbl.size() != nInp)
        {
          Lucee::Except lce("PointSourceIfc::readInput: The table 'inpComponents' must have ");
          lce << nInp << " entries. Instead provided " << inpListDbl.size();
          throw lce;
        }
        else
        {
          for (unsigned i=0; i<nInp; ++i)
            inpComponents[i] = (unsigned) inpListDbl[i];
        }
      }
      else
      {
// resize and set input variable list from one specified in Lua
        nInp = inpListDbl.size();
        inpComponents.resize(nInp);
        for (unsigned i=0; i<nInp; ++i)
          inpComponents[i] = (unsigned) inpListDbl[i];
      }
    }

    if (tbl.hasNumVec("outComponents"))
    {
// get output component list
      std::vector<double> outListDbl = tbl.getNumVec("outComponents");
      if (allowArb == false)
      {
        if (outListDbl.size() != nOut)
        {
          Lucee::Except lce("PointSourceIfc::readInput: The table 'outComponents' must have ");
          lce << nOut << " entries. Instead provided " << outListDbl.size();
          throw lce;
        }
        else
        {
          for (unsigned i=0; i<nOut; ++i)
            outComponents[i] = (unsigned) outListDbl[i];
        }
      }
      else
      {
// resize and set output variable list from one specified in Lua
        nOut = outListDbl.size();
        outComponents.resize(nOut);
        for (unsigned i=0; i<nOut; ++i)
          outComponents[i] = (unsigned) outListDbl[i];
      }
    }
// allocate space needed for computing source
    out.resize(nOut);
// allocate space needed for computing source jacobian
    srcJac = Lucee::Matrix<double>(nOut, nOut);
  }

  void
  PointSourceIfc::calcSource(double tm, const double loc[3], const double *inp, double *src)
  {
// set data pointer so derived classes can get needed variables
    data = inp;
    for (unsigned i=0; i<nOut; ++i) out[i] = 0.0;
// call derived class method to compute source
    this->getSource(tm, loc, out);
// copy source over into proper location
    for (unsigned i=0; i<nOut; ++i)
      src[outComponents[i]] = out[i];
  }

  void
  PointSourceIfc::calcSourceJac(double tm, const double loc[3], const double *inp, 
    Lucee::Matrix<double>& jac)
  {
// set data pointer so derived classes can get needed variables
    data = inp;
// clear contents of jacobian before getting it from derived class
    srcJac = 0.0;
// call derived class method to compute jacobian
    this->getSourceJac(tm, loc, srcJac);
// copy over into proper location
    for (unsigned i=0; i<nOut; ++i)
      for (unsigned j=0; j<nOut; ++j)
        jac(outComponents[i], outComponents[j]) = srcJac(i,j);
  }

  void
  PointSourceIfc::getSourceJac(double tm, const double loc[3], 
    Lucee::Matrix<double>& jac)
  {
    throw Lucee::Except("PointSourceIfc::getSourceJac: Not implemented!");
  }
}
