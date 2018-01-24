/**
 * @file	LcGeneralZonalAverage.cpp
 *
 * @brief       Updater to compute zonal average of a ndim field/distf.
 * Output will be a ndim-1 field.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGeneralZonalAverage.h>
#include <LcGlobals.h>
#include <LcLinAlgebra.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <cmath>
#include <vector>

namespace Lucee
{
  //template <> const char *GeneralZonalAverage<1>::id = "GeneralZonalAverage1D";
  template <> const char *GeneralZonalAverage<2>::id = "GeneralZonalAverage2D";
  template <> const char *GeneralZonalAverage<3>::id = "GeneralZonalAverage3D";
  template <> const char *GeneralZonalAverage<4>::id = "GeneralZonalAverage4D";
  template <> const char *GeneralZonalAverage<5>::id = "GeneralZonalAverage5D";

  template <unsigned NDIM>
  bool
  GeneralZonalAverage<NDIM>::sameAveCoords(unsigned n, unsigned cn, double dxMin,
    unsigned srcDir[NDIM-1], const Eigen::MatrixXd& phaseC, const Eigen::MatrixXd& confC)
  {
    const unsigned ADIM = NDIM-1;
    for (unsigned dim = 0; dim<ADIM; ++dim)
      if (! (std::fabs(phaseC(n, srcDir[dim])-confC(cn, dim))<1e-4*dxMin) )
        return false;
    return true;
  }

  template <unsigned NDIM>
  GeneralZonalAverage<NDIM>::GeneralZonalAverage()
    : Lucee::UpdaterIfc()
  {
    momentLocal = 0;
  }

  template <unsigned NDIM>
  GeneralZonalAverage<NDIM>::~GeneralZonalAverage()
  {
    delete momentLocal;
  }

  template <unsigned NDIM>
  void
  GeneralZonalAverage<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);
    const unsigned ADIM = NDIM-1;

    // get hold of element to use for source field
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("sourceBasis"))
      sourceBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("sourceBasis");
    else
      throw Lucee::Except(
        "GeneralZonalAverage::readInput: Must specify element to use using 'sourceBasis'");

    // get hold of element to use for average field
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<ADIM> >("aveBasis"))
      aveBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<ADIM> >("aveBasis");
    else
      throw Lucee::Except(
        "GeneralZonalAverage::readInput: Must specify element to use using 'aveBasis'");

    intDir = 1;
  }

  template <unsigned NDIM>
  void
  GeneralZonalAverage<NDIM>::initialize()
  {
    Lucee::UpdaterIfc::initialize();
    const unsigned ADIM = NDIM-1;

    std::vector<int> coordinateMap;
    for (int i = 0; i < ADIM; i++)
        coordinateMap.push_back(i);
    
    // get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    // local region to update
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    // // Add something like this? (from LcCopyNodalFields)
    // Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    // seq.step(); // just to get to first index
    // int idx[NDIM];
    // seq.fillWithIndex(idx);
    // sourceBasis->setIndex(idx);
    // int idxAve[ADIM];
    // for (int aIndex = 0; aIndex < ADIM; aIndex++)
    //   idxAve[aIndex] = idx[coordinateMap[aIndex]];
    // sourceBasis->setIndex(idxAve); // only first DIM elements are used
    
    // get number of nodes in initial field and averaged field
    unsigned nlocalAve = aveBasis->getNumNodes();
    unsigned nlocalSource = sourceBasis->getNumNodes();

    // get volume interpolation matrices for averaged field element
    int nVolQuadAve = aveBasis->getNumGaussNodes();
    std::vector<double> volWeightsAve(nVolQuadAve);
    Lucee::Matrix<double> tempVolQuadAve(nVolQuadAve, nlocalAve);
    Lucee::Matrix<double> tempVolCoordsAve(nVolQuadAve, ANC);

    aveBasis->getGaussQuadData(tempVolQuadAve, tempVolCoordsAve, volWeightsAve);

    Eigen::MatrixXd volQuadAve(nVolQuadAve, nlocalAve);
    copyLuceeToEigen(tempVolQuadAve, volQuadAve);
    Eigen::MatrixXd volCoordsAve(nVolQuadAve, (unsigned) ANC);
    copyLuceeToEigen(tempVolCoordsAve, volCoordsAve);

    
    // get volume interpolation matrices for source field element
    int nVolQuadSource = sourceBasis->getNumGaussNodes();
    std::vector<double> volWeightsSource(nVolQuadSource);
    Lucee::Matrix<double> tempVolQuadSource(nVolQuadSource, nlocalSource);
    Lucee::Matrix<double> tempVolCoordsSource(nVolQuadSource, (unsigned) SNC);

    sourceBasis->getGaussQuadData(tempVolQuadSource, tempVolCoordsSource, volWeightsSource);

    Eigen::MatrixXd volQuadSource(nVolQuadSource, nlocalSource);
    copyLuceeToEigen(tempVolQuadSource, volQuadSource);
    Eigen::MatrixXd volCoordsSource(nVolQuadSource, (unsigned) SNC);
    copyLuceeToEigen(tempVolCoordsSource, volCoordsSource);

    // compute mapping of source field quadrature nodes to averaged field space
    // quadrature nodes. The assumption here is that the node layout in source space
    // and ave space are such that each node in source space has
    // exactly one node co-located with it in ave space. No
    // "orphan" source space nodes are allowed, and an exception is thrown
    // if that occurs.
    srcAveMap.resize(volWeightsSource.size());
    unsigned srcDir[ADIM];
    unsigned count = 0;
    for (unsigned d=0; d<NDIM; ++d)
    {
      if (d != intDir) 
        srcDir[count++] = d;
    }
    
    double dxMin = grid.getDx(0);
    for (unsigned d=1; d<ADIM; ++d)
      dxMin = std::min(dxMin, grid.getDx(d));

    for (unsigned n=0; n<volWeightsSource.size(); ++n)
    {
      bool saFound = false;
      for (unsigned cn=0; cn<volWeightsAve.size(); ++cn)
        if (sameAveCoords(n, cn, dxMin, srcDir, volCoordsSource, volCoordsAve))
        {
          srcAveMap[n] = cn;
          saFound = true;
          break;
        }
      if (!saFound)
      {
        Lucee::Except lce(
          "GeneralZonalAverage::initialize: No matching averaged-field quadrature node for source-field quadrature node ");
        lce << n;
        throw lce;
      }
    }

    // Get Averaged Field Mass Matrix
    Lucee::Matrix<double> tempMassMatrixAve(nlocalAve, nlocalAve);
    aveBasis->getMassMatrix(tempMassMatrixAve);
    Eigen::MatrixXd massMatrixAve(nlocalAve, nlocalAve);
    copyLuceeToEigen(tempMassMatrixAve, massMatrixAve);

    mom0Matrix = Eigen::MatrixXd::Zero(nlocalAve, nlocalSource);

    for (int i=0; i<nlocalAve; ++i)
    {
      for (int j=0; j<nlocalSource; ++j)
      {
        double integralResult = 0.0;
        // Compute integral of phiAve_i * phi_j
        for (int gaussIndex = 0; gaussIndex < volWeightsSource.size(); ++gaussIndex)
          integralResult
            += volWeightsSource[gaussIndex]*volQuadAve(srcAveMap[gaussIndex],i)*volQuadSource(gaussIndex,j);
        mom0Matrix(i, j) = integralResult;
      }
    }

    // Multiply matrices by inverse of mass matrix
    mom0Matrix = massMatrixAve.inverse()*mom0Matrix;
    //std::cout << " mom0Matrix " << mom0Matrix << std::endl;
    
  }
   
  template <unsigned NDIM>
  Lucee::UpdaterStatus
  GeneralZonalAverage<NDIM>::update(double t)
  {
    const unsigned ADIM = NDIM-1;
    
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    // get input field (NDIM)
    const Lucee::Field<NDIM, double>& fld = this->getInp<Lucee::Field<NDIM, double> >(0);
    //std::cout << " input field " << fld << std::endl;
    
    // get output field (ADIM = NDIM-1)
    Lucee::Field<ADIM, double>& momentGlobal = this->getOut<Lucee::Field<ADIM, double> >(0);
    
    std::vector<int> coordinateMap;
    for (int i = 0; i < ADIM; i++)
      coordinateMap.push_back(i);

    if (!momentLocal)
    {
      // allocate memory for local moment calculation if not already
      // done: we need to ensure space is also allocated for the
      // ghost-cells as otherwise there is a size mis-match in the
      // allReduce call to sync across velocity space
      Lucee::Region<ADIM, int> localRgn = momentGlobal.getRegion();
      Lucee::Region<ADIM, int> localExtRgn = momentGlobal.getExtRegion();
      
      int lowerAve[ADIM], upperAve[ADIM], lg[ADIM], ug[ADIM];
      for (int i=0; i<ADIM; ++i)
      {
        lowerAve[i] = localRgn.getLower(i);
        upperAve[i] = localRgn.getUpper(i);
        lg[i] = localRgn.getLower(i) - localExtRgn.getLower(i);
        ug[i] = localExtRgn.getUpper(i) - localRgn.getUpper(i);
      }
      Lucee::Region<ADIM, int> rgnAve(lowerAve, upperAve);
      momentLocal = new Lucee::Field<ADIM, double>(rgnAve, momentGlobal.getNumComponents(), lg, ug);
    }

    // local phase-space region to update
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();

    int idx[NDIM];
    int idxAve[ADIM];
    double xc[SNC];
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    unsigned nlocalAve = aveBasis->getNumNodes();
    unsigned nlocalSource = sourceBasis->getNumNodes();

    // total number of local averaged-field space cells
    int localPositionCells = momentGlobal.getExtRegion().getVolume();

    // clear out contents of output field
    (*momentLocal) = 0.0;

    unsigned srcDir[ADIM];
    unsigned count = 0;
    for (unsigned d=0; d<NDIM; ++d)
    {
      if (d != intDir) 
        srcDir[count++] = d;
    }    

    // iterators into fields
    Lucee::ConstFieldPtr<double> fldPtr = fld.createConstPtr();
    Lucee::FieldPtr<double> momentPtr = momentLocal->createPtr();
    Eigen::VectorXd fldVec(nlocalSource);
     
    while(seq.step())
    {
      seq.fillWithIndex(idx);
      for (unsigned d=0; d<ADIM; ++d)
        idxAve[d] = idx[srcDir[d]];
      
      grid.setIndex(idx);
      fld.setPtr(fldPtr, idx);
      momentLocal->setPtr(momentPtr, idxAve);
      
      for (int i = 0; i < nlocalSource; i++)
        fldVec(i) = fldPtr[i];

      //std::cout << "fldVec in seq" << fldVec << std::endl;
      // Accumulate contribution to moment from this cell
      Eigen::VectorXd resultVector = mom0Matrix*fldVec;

      for (int i = 0; i < nlocalAve; i++)
        momentPtr[i] += resultVector(i);
    }

    // we need to get moment communicator of field as updater's moment
    // communicator is same as its grid's moment communicator. In this
    // case, grid is phase-space grid, which is not what we want.
    TxCommBase *momComm = momentGlobal.getMomComm();
    unsigned xsize = localPositionCells*nlocalAve; // amount to communicate, cells*nodes_per_cell
    momComm->allreduce(xsize, &momentLocal->first(), &momentGlobal.first(), TX_SUM);
  
    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  GeneralZonalAverage<NDIM>::declareTypes()
  {
    const unsigned ADIM = NDIM-1;
    // Input potential on a NDIM field
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    // Zonal-averaged potential on a NDIM-1 field
    this->appendOutVarType(typeid(Lucee::Field<ADIM, double>));
  }

  template <unsigned NDIM>
  void
  GeneralZonalAverage<NDIM>::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }

  //instantiations
  //template class GeneralZonalAverage<1>;
  template class GeneralZonalAverage<2>;
  template class GeneralZonalAverage<3>;
  template class GeneralZonalAverage<4>;
  template class GeneralZonalAverage<5>;

}
