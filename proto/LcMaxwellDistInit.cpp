/**
 * @file	LcMaxwellDistInit.cpp
 *
 * @brief	Updater for initializin the Maxwellian distribution from moments
 */

// lucee includes
#include <LcMaxwellDistInit.h>

namespace Lucee
{
// set ids for module system
  template <> const char *MaxwellDistInit<1, 1>::id = "MaxwellDistInit1X1V";
  template <> const char *MaxwellDistInit<1, 2>::id = "MaxwellDistInit1X2V";
  template <> const char *MaxwellDistInit<1, 3>::id = "MaxwellDistInit1X3V";
  template <> const char *MaxwellDistInit<2, 2>::id = "MaxwellDistInit2X2V";
  template <> const char *MaxwellDistInit<2, 3>::id = "MaxwellDistInit2X3V";
  //template <> const char *MaxwellDistInit<3, 3>::id = "MaxwellDistInit3X3V";

  template <unsigned CDIM, unsigned VDIM>
  bool
  MaxwellDistInit<CDIM,VDIM>::sameConfigCoords(unsigned n, unsigned cn, double dxMin,
    const Lucee::Matrix<double>& phaseC, const Lucee::Matrix<double>& confC)
  {
    for (unsigned dim = 0; dim<CDIM; ++dim)
      if (! (std::fabs(phaseC(n, dim)-confC(cn, dim))<1e-4*dxMin) )
        return false;
    return true;
  }
  //----------------------------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  MaxwellDistInit<CDIM,VDIM>::MaxwellDistInit()
    : UpdaterIfc()
  {
  }
  template <unsigned CDIM, unsigned VDIM>
  MaxwellDistInit<CDIM, VDIM>::~MaxwellDistInit()
  {
  }
  //----------------------------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  void
  MaxwellDistInit<CDIM, VDIM>::readInput(Lucee::LuaTable& tbl)
  {
    const unsigned NDIM = CDIM+VDIM;
    // call base class method
    UpdaterIfc::readInput(tbl);

    // get hold on the basis
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("phaseBasis"))
      phaseBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("phaseBasis");
    else
      throw Lucee::Except("MaxwellDistInit::readInput: Must specify phase space basis to use using 'phaseBasis'");
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<CDIM> >("confBasis"))
      confBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<CDIM> >("confBasis");
    else
      throw Lucee::Except("MaxwellDistInit::readInput: Must specify configuration space basis to use using 'confBasis'");

    // flag for outputing Maxwellian with zero drift velocity
    zeroDriftOutput = false;
    if (tbl.hasBool("zeroDriftOutput"))
      zeroDriftOutput = tbl.getBool("zeroDriftOutput");

    // flag to output Maxwellian with arbitrary density (original
    // number density is still need for calculation of u and vt from
    // first and second moment
    arbitraryDensity = false;
    if (tbl.hasBool("arbitraryDensity"))
      arbitraryDensity = tbl.getBool("arbitraryDensity");
  }
  //----------------------------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  void 
  MaxwellDistInit<CDIM, VDIM>::initialize()
  {
    UpdaterIfc::initialize();
 
    const unsigned NDIM = CDIM+VDIM;
    unsigned nlocal = phaseBasis->getNumNodes();

    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[NDIM];
    seq.fillWithIndex(idx);
    phaseBasis->setIndex(idx);
    confBasis->setIndex(idx); // only first CDIM elements are used

    // compute mapping of phase-space nodes to configuration space
    // nodes. The assumption here is that the node layout in phase-space
    // and configuration space are such that each node in phase-space has
    // exactly one node co-located with it in configuration space. No
    // "orphan" phase-space node are allowed, and an exception is thrown
    // if that occurs.
    phaseConfMap.resize(nlocal);
    Lucee::Matrix<double> phaseNodeCoords(phaseBasis->getNumNodes(), PNC);
    Lucee::Matrix<double> confNodeCoords(confBasis->getNumNodes(), CNC);

    double dxMin = grid.getDx(0);
    for (unsigned d=1; d<CDIM; ++d)
      dxMin = std::min(dxMin, grid.getDx(d));

    phaseBasis->getNodalCoordinates(phaseNodeCoords);
    confBasis->getNodalCoordinates(confNodeCoords);
    for (unsigned n=0; n<nlocal; ++n)
    {
      bool pcFound = false;
      for (unsigned cn=0; cn<confBasis->getNumNodes(); ++cn)
        if (sameConfigCoords(n, cn, dxMin, phaseNodeCoords, confNodeCoords))
        {
          phaseConfMap[n] = cn;
          pcFound = true;
          break;
        }
      if (!pcFound)
      {
        Lucee::Except lce(
          "MaxwellDistInit::initialize: No matching configuration space node for phase-space node ");
        lce << n;
        throw lce;
      }
    }
  }
  //----------------------------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  Lucee::UpdaterStatus 
  MaxwellDistInit<CDIM, VDIM>::update(double t)
  {
    const unsigned NDIM = CDIM+VDIM;

    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    const Lucee::Field<CDIM, double>& zerothMoment =
      this->getInp<Lucee::Field<CDIM, double> >(0);
    const Lucee::Field<CDIM, double>& firstMoment =
      this->getInp<Lucee::Field<CDIM, double> >(1);
    const Lucee::Field<CDIM, double>& secondMoment =
      this->getInp<Lucee::Field<CDIM, double> >(2);
      
    // For compatibility with mathematical definition of DG, I am
    // retaining the names "q" for the distribution function.
    Lucee::Field<NDIM, double>& q =
      this->getOut<Lucee::Field<NDIM, double> >(0);

    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    Lucee::Region<NDIM, double> compSpace = grid.getComputationalSpace();

    Lucee::FieldPtr<double> qPtr = q.createPtr();
    Lucee::ConstFieldPtr<double> zerothMomentPtr =
      zerothMoment.createConstPtr();
    Lucee::ConstFieldPtr<double> firstMomentPtr =
      firstMoment.createConstPtr();
    Lucee::ConstFieldPtr<double> secondMomentPtr =
      secondMoment.createConstPtr();

    Lucee::Matrix<double> phaseNodeCoords(phaseBasis->getNumNodes(), PNC);

    q = 0.0; // use q to store increment initially    
    
    int idx[NDIM];
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);

    unsigned numNodesPhase = phaseBasis->getNumNodes();
    unsigned numNodesConf = confBasis->getNumNodes();
    
    std::vector<double> dens(numNodesConf);
    std::vector<double> densInv(numNodesConf);
    double densIn;
    std::vector<double> vTerm2(numNodesConf);
    double vTerm2In;
    Lucee::Matrix<double> vDrift(numNodesConf, VDIM);
    double v[VDIM];

    // get an "area" of velocity space so the updater can return a
    // uniform distribution with correct density in case the higher
    // moments were not calculated correctly and vTerm2 < 0
    double vSpaceArea = 1;
    double vSpaceAreaInv;
    for (int dim = 0; dim<VDIM; ++dim)
    {
      vSpaceArea *= compSpace.getShape(CDIM+dim);
    }
    vSpaceAreaInv = 1/vSpaceArea;

    if (!arbitraryDensity)
    {
      while (seq.step())
      {
	seq.fillWithIndex(idx);
	q.setPtr(qPtr, idx);
	zerothMoment.setPtr(zerothMomentPtr, idx);
	firstMoment.setPtr(firstMomentPtr, idx);
	secondMoment.setPtr(secondMomentPtr, idx);

	phaseBasis->setIndex(idx);
	phaseBasis->getNodalCoordinates(phaseNodeCoords);
      
	for (unsigned nodeIdx = 0; nodeIdx<numNodesConf; ++nodeIdx)
	{
	  dens[nodeIdx] = zerothMomentPtr[nodeIdx];
	  densInv[nodeIdx] = 1/dens[nodeIdx];
	  vTerm2[nodeIdx] = densInv[nodeIdx]*secondMomentPtr[nodeIdx];
	  for (unsigned dim = 0; dim<VDIM; ++dim)
	  {
	    vDrift(nodeIdx, dim) = densInv[nodeIdx]*
	      firstMomentPtr[nodeIdx*VDIM+dim];
	    vTerm2[nodeIdx] = vTerm2[nodeIdx]-
	      vDrift(nodeIdx, dim)*vDrift(nodeIdx, dim);
	  }  	
	}
	for (unsigned nodeIdx = 0; nodeIdx<numNodesPhase; ++nodeIdx) 
	{
	  if (zeroDriftOutput) 
	    for (unsigned dim = 0; dim<VDIM; ++dim)
	      v[dim] = phaseNodeCoords(nodeIdx, CDIM+dim);
	  else
	    for (unsigned dim = 0; dim<VDIM; ++dim)
	      v[dim] = phaseNodeCoords(nodeIdx, CDIM+dim)-
		vDrift(phaseConfMap[nodeIdx], dim);
	  densIn = dens[phaseConfMap[nodeIdx]];

	  vTerm2In = vTerm2[phaseConfMap[nodeIdx]];
	  if (vTerm2In >= 0) 
	    qPtr[nodeIdx] = evaluateMaxwell(v, densIn, vTerm2In);
	  else
	    qPtr[nodeIdx] = densIn*vSpaceAreaInv;
	}
      }
    }
    else
    {
      // declaration and usage of this field must be in one block
      const Lucee::Field<CDIM, double>& zerothMomentOut =
	this->getInp<Lucee::Field<CDIM, double> >(3);
      Lucee::ConstFieldPtr<double> zerothMomentOutPtr = 
	zerothMomentOut.createConstPtr();

      while (seq.step())
      {
	seq.fillWithIndex(idx);
	q.setPtr(qPtr, idx);
	zerothMoment.setPtr(zerothMomentPtr, idx);
	firstMoment.setPtr(firstMomentPtr, idx);
	secondMoment.setPtr(secondMomentPtr, idx);
	zerothMomentOut.setPtr(zerothMomentOutPtr, idx);

	phaseBasis->setIndex(idx);
	phaseBasis->getNodalCoordinates(phaseNodeCoords);
      
	for (unsigned nodeIdx = 0; nodeIdx<numNodesConf; ++nodeIdx)
	{
	  dens[nodeIdx] = zerothMomentPtr[nodeIdx];
	  densInv[nodeIdx] = 1/dens[nodeIdx];
	  vTerm2[nodeIdx] = densInv[nodeIdx]*secondMomentPtr[nodeIdx];
	  for (unsigned dim = 0; dim<VDIM; ++dim)
	  {
	    vDrift(nodeIdx, dim) = densInv[nodeIdx]*
	      firstMomentPtr[nodeIdx*VDIM+dim];
	    vTerm2[nodeIdx] = vTerm2[nodeIdx]-
	      vDrift(nodeIdx, dim)*vDrift(nodeIdx, dim);
	  }  	
	}
	for (unsigned nodeIdx = 0; nodeIdx<numNodesPhase; ++nodeIdx) 
	{
	  if (zeroDriftOutput) 
	    for (unsigned dim = 0; dim<VDIM; ++dim)
	      v[dim] = phaseNodeCoords(nodeIdx, CDIM+dim);
	  else
	    for (unsigned dim = 0; dim<VDIM; ++dim)
	      v[dim] = phaseNodeCoords(nodeIdx, CDIM+dim)-
		vDrift(phaseConfMap[nodeIdx], dim);
	  densIn = 1/zerothMomentOutPtr[phaseConfMap[nodeIdx]];

	  vTerm2In = vTerm2[phaseConfMap[nodeIdx]];
	  if (vTerm2In >= 0) 
	    qPtr[nodeIdx] = evaluateMaxwell(v, densIn, vTerm2In);
	  else
	    qPtr[nodeIdx] = densIn*vSpaceAreaInv;
	}
      }
    }

    return Lucee::UpdaterStatus();
  }
  //----------------------------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  void MaxwellDistInit<CDIM, VDIM>::declareTypes()
  {
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }
  //----------------------------------------------------------------------------
  template <unsigned CDIM, unsigned VDIM>
  double MaxwellDistInit<CDIM, VDIM>::evaluateMaxwell(double v[VDIM], double n, 
    double vTerm2)
  {
    double result = 0;
    double factor = 1/sqrt(2*M_PI*vTerm2);
    for (unsigned dim = 0; dim<VDIM; ++dim)
      result += v[dim]*v[dim];
    result = n*exp(-0.5*result/vTerm2);
    for (unsigned dim = 0; dim<VDIM; ++dim)
      result *= factor;
    return result;
  }
  //----------------------------------------------------------------------------
  // instantiations
  template class MaxwellDistInit<1, 1>;
  template class MaxwellDistInit<1, 2>;
  template class MaxwellDistInit<1, 3>;
  template class MaxwellDistInit<2, 2>;
  template class MaxwellDistInit<2, 3>;
  //template class MaxwellDistInit<3, 3>;
}

