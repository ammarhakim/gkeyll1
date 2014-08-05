/**
 * @file	LcCompletePolynomialElement.cpp
 *
 * @brief Complete polynomial element
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcExcept.h>
#include <LcGridIfc.h>
#include <LcCompletePolynomialElement.h>

namespace Lucee
{
  using namespace Eigen;

// set module name
  template <> const char *CompletePolynomialElement<1>::id = "CompletePolynomialElement";
  template <> const char *CompletePolynomialElement<2>::id = "CompletePolynomialElement";
  template <> const char *CompletePolynomialElement<3>::id = "CompletePolynomialElement";
  template <> const char *CompletePolynomialElement<4>::id = "CompletePolynomialElement";
  template <> const char *CompletePolynomialElement<5>::id = "CompletePolynomialElement";

  template <unsigned NDIM>
  CompletePolynomialElement<NDIM>::CompletePolynomialElement()
    : NodalFiniteElementIfc<NDIM>(1)
  {

  }
  
  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::NodalFiniteElementIfc<NDIM>::readInput(tbl);
    polyOrder = (unsigned) tbl.getNumber("polyOrder");
    // Compute maximum polynomial power we want (no transformations)
    // This will be the maximum power of a single variable in anything
    // evaluated (e.g. product of basis functions)
    // Note: don't change this without paying attention to how maxPower is
    // used in the rest of the code! Affects more than just gaussPoints.
    if (polyOrder == 1)
      maxPower = 4;
    else
      maxPower = 3*polyOrder;

    this->setNumNodes(getCompletePolynomialDimension(polyOrder, NDIM));
    setupMatrices();
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::getExclusiveNodeIndices(std::vector<int>& ndIds)
  {
    throw Lucee::Except("CompletePolynomialElement::getExclusiveNodeIndices: not implemented.");
  }

  template <unsigned NDIM>
  unsigned
  CompletePolynomialElement<NDIM>::getNumSurfLowerNodes(unsigned dir) const
  {
    return this->getNumNodes();
  }

  template <unsigned NDIM>
  unsigned
  CompletePolynomialElement<NDIM>::getNumSurfUpperNodes(unsigned dir) const
  {
    return this->getNumNodes();
  }

  template <unsigned NDIM>
  unsigned
  CompletePolynomialElement<NDIM>::getNumGlobalNodes() const
  {
    throw Lucee::Except("CompletePolynomialElement::getNumGlobalNodes: not implemented.");
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::getLocalToGlobal(std::vector<int>& lgMap) const
  {
    throw Lucee::Except("CompletePolynomialElement::getLocalToGlobal: not implemented.");
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::getSurfLowerLocalToGlobal(unsigned dim,
    std::vector<int>& lgMap) const
  {
    throw Lucee::Except("CompletePolynomialElement::getSurfLowerLocalToGlobal: not implemented.");
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::getSurfUpperLocalToGlobal(unsigned dim,
    std::vector<int>& lgMap) const
  {
    throw Lucee::Except("CompletePolynomialElement::getSurfUpperLocalToGlobal: not implemented.");
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::getSurfLowerNodeNums(unsigned dir,
    std::vector<int>& nodeNum) const
  {
    for (int i = 0; i < this->getNumNodes(); i++)
      nodeNum[i] = i;
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::getSurfUpperNodeNums(unsigned dir,
    std::vector<int>& nodeNum) const
  {
    for (int i = 0; i < this->getNumNodes(); i++)
      nodeNum[i] = i;
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::getNodalCoordinates(Lucee::Matrix<double>& nodeCoords)
  {
    throw Lucee::Except("CompletePolynomialElement::getNodalCoordinates: not implemented.");
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::getWeights(std::vector<double>& w)
  {
    throw Lucee::Except("CompletePolynomialElement::getWeights: Not implemented!");
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::getSurfUpperWeights(unsigned dir, std::vector<double>& w)
  {
    throw Lucee::Except("CompletePolynomialElement::getSurfUpperWeights: Not implemented!");
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::getSurfLowerWeights(unsigned dir, std::vector<double>& w)
  {
    throw Lucee::Except("CompletePolynomialElement::getSurfLowerWeights: Not implemented!");
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::getMassMatrix(Lucee::Matrix<double>& NjNk) const
  {
    for (int i = 0; i < refMass.rows(); i++)
      for (int j = 0; j < refMass.cols(); j++)
        NjNk(i, j) = refMass(i, j);
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::getLowerFaceMassMatrix(unsigned dir,
    Lucee::Matrix<double>& NjNk) const
  {
    for (int i = 0; i < refFaceMassLower[dir].rows(); i++)
      for (int j = 0; j < refFaceMassLower[dir].cols(); j++)
        NjNk(i, j) = refFaceMassLower[dir](i, j);
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::getUpperFaceMassMatrix(unsigned dir,
    Lucee::Matrix<double>& NjNk) const
  {
    for (int i = 0; i < refFaceMassUpper[dir].rows(); i++)
      for (int j = 0; j < refFaceMassUpper[dir].cols(); j++)
        NjNk(i, j) = refFaceMassUpper[dir](i, j);
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::getStiffnessMatrix(Lucee::Matrix<double>& DNjDNk) const
  {
    for (int i = 0; i < refStiffness.rows(); i++)
      for (int j = 0; j < refStiffness.cols(); j++)
        DNjDNk(i, j) = refStiffness(i, j);
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::getGradStiffnessMatrix(
    unsigned dir, Lucee::Matrix<double>& DNjNk) const
  {
    for (int i = 0; i < refGradStiffness[dir].rows(); i++)
      for (int j = 0; j < refGradStiffness[dir].cols(); j++)
        DNjNk(i, j) = refGradStiffness[dir](i, j);
  }

  template <unsigned NDIM>
  unsigned
  CompletePolynomialElement<NDIM>::getNumGaussNodes() const
  {
    return gaussNodeList.rows();
  }

  template <unsigned NDIM>
  unsigned
  CompletePolynomialElement<NDIM>::getNumSurfGaussNodes() const
  {
    return gaussNodeListLowerSurf[0].rows();
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::getGaussQuadData(Lucee::Matrix<double>& interpMat,
    Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const
  {
    double weightScale = 1.0;

    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
      weightScale *= 0.5*dq[dimIndex];

    for (int gaussIndex = 0; gaussIndex < gaussNodeList.rows(); gaussIndex++)
    {
      if (NDIM == 1)
      {
        ordinates(gaussIndex,0) = gaussNodeList(gaussIndex,0);
        ordinates(gaussIndex,1) = 0;
        ordinates(gaussIndex,2) = 0;
      }
      else if (NDIM == 2)
      {
        ordinates(gaussIndex,0) = gaussNodeList(gaussIndex,0);
        ordinates(gaussIndex,1) = gaussNodeList(gaussIndex,1);
        ordinates(gaussIndex,2) = 0;
      }
      else
      {
        for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
          ordinates(gaussIndex,dimIndex) = gaussNodeList(gaussIndex,dimIndex);
      }

      weights[gaussIndex] = gaussNodeList(gaussIndex,NDIM)*weightScale;

      // Copy interpolation matrix
      for (int basisIndex = 0; basisIndex < functionEvaluations.cols(); basisIndex++)
        interpMat(gaussIndex, basisIndex) = functionEvaluations(gaussIndex, basisIndex);
    }
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::getSurfLowerGaussQuadData(unsigned dir, Lucee::Matrix<double>& interpMat,
    Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const
  {
    double weightScale = 1.0;

    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      if (dimIndex != dir)
        weightScale *= 0.5*dq[dimIndex];
    }

    for (int gaussIndex = 0; gaussIndex < gaussNodeListLowerSurf[dir].rows(); gaussIndex++)
    {
      if (NDIM == 1)
      {
        ordinates(gaussIndex,0) = gaussNodeListLowerSurf[dir](gaussIndex,0);
        ordinates(gaussIndex,1) = 0;
        ordinates(gaussIndex,2) = 0;
      }
      else if (NDIM == 2)
      {
        ordinates(gaussIndex,0) = gaussNodeListLowerSurf[dir](gaussIndex,0);
        ordinates(gaussIndex,1) = gaussNodeListLowerSurf[dir](gaussIndex,1);
        ordinates(gaussIndex,2) = 0;
      }
      else
      {
        for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
          ordinates(gaussIndex,dimIndex) = gaussNodeListLowerSurf[dir](gaussIndex,dimIndex);
      }

      weights[gaussIndex] = gaussNodeListLowerSurf[dir](gaussIndex,NDIM)*weightScale;

      // Copy interpolation matrix
      for (int basisIndex = 0; basisIndex < lowerSurfaceEvaluations[dir].cols(); basisIndex++)
        interpMat(gaussIndex, basisIndex) = lowerSurfaceEvaluations[dir](gaussIndex, basisIndex);
    }
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::getSurfUpperGaussQuadData(unsigned dir, Lucee::Matrix<double>& interpMat,
    Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const
  {
    double weightScale = 1.0;

    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      if (dimIndex != dir)
        weightScale *= 0.5*dq[dimIndex];
    }

    for (int gaussIndex = 0; gaussIndex < gaussNodeListUpperSurf[dir].rows(); gaussIndex++)
    {
      if (NDIM == 1)
      {
        ordinates(gaussIndex,0) = gaussNodeListUpperSurf[dir](gaussIndex,0);
        ordinates(gaussIndex,1) = 0;
        ordinates(gaussIndex,2) = 0;
      }
      else if (NDIM == 2)
      {
        ordinates(gaussIndex,0) = gaussNodeListUpperSurf[dir](gaussIndex,0);
        ordinates(gaussIndex,1) = gaussNodeListUpperSurf[dir](gaussIndex,1);
        ordinates(gaussIndex,2) = 0;
      }
      else
      {
        for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
          ordinates(gaussIndex,dimIndex) = gaussNodeListUpperSurf[dir](gaussIndex,dimIndex);
      }

      weights[gaussIndex] = gaussNodeListUpperSurf[dir](gaussIndex,NDIM)*weightScale;

      // Copy interpolation matrix
      for (int basisIndex = 0; basisIndex < upperSurfaceEvaluations[dir].cols(); basisIndex++)
        interpMat(gaussIndex, basisIndex) = upperSurfaceEvaluations[dir](gaussIndex, basisIndex);
    }
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::getMomentMatrix(unsigned p, Lucee::Matrix<double>& momMat) const
  {
    if (p > maxMoment || NDIM != 2)
    {
      // moments higher than 3 not supported for now
      Lucee::Except lce("CompletePolynomialElement::getMomentMatrix: Moment matrix of order ");
      lce << p << " not supported. NDIM must also be 2.";
      throw lce;
    }
    else
    {
      for (int rowIndex = 0; rowIndex < momMatrix[p].rows(); rowIndex++)
        for (int colIndex = 0; colIndex < momMatrix[p].cols(); colIndex++)
          momMat(rowIndex,colIndex) = momMatrix[p](rowIndex,colIndex);
    }
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::getDiffusionMatrices(std::vector<Lucee::Matrix<double> >& iMat_o,
    std::vector<Lucee::Matrix<double> >& lowerMat_o, std::vector<Lucee::Matrix<double> >& upperMat_o) const
  {
    if (polyOrder > 2)
      Lucee::Except("CompletePolynomialElement::getDiffusionMatrices: Not implemented for higher than quadratic!");
    if (NDIM != 2)
      Lucee::Except("CompletePolynomialElement::getDiffusionMatrices: Only implemented for 2D");

    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      for (int rowIndex = 0; rowIndex < iMatDiffusion[dimIndex].rows(); rowIndex++)
      {
        for (int colIndex = 0; colIndex < iMatDiffusion[dimIndex].cols(); colIndex++)
        {
          iMat_o[dimIndex](rowIndex,colIndex) = iMatDiffusion[dimIndex](rowIndex,colIndex);
          lowerMat_o[dimIndex](rowIndex,colIndex) = lowerMatDiffusion[dimIndex](rowIndex,colIndex);
          upperMat_o[dimIndex](rowIndex,colIndex) = upperMatDiffusion[dimIndex](rowIndex,colIndex);
        }
      }
    }
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::getHyperDiffusionMatrices(std::vector<Lucee::Matrix<double> >& iMat_o,
    std::vector<Lucee::Matrix<double> >& lowerMat_o, std::vector<Lucee::Matrix<double> >& upperMat_o) const
  {
    if (polyOrder > 2)
      Lucee::Except("CompletePolynomialElement::getHyperDiffusionMatrices: Not implemented for higher than quadratic!");
    if (NDIM != 2)
      Lucee::Except("CompletePolynomialElement::getHyperDiffusionMatrices: Only implemented for 2D");

    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      for (int rowIndex = 0; rowIndex < iMatHyperDiffusion[dimIndex].rows(); rowIndex++)
      {
        for (int colIndex = 0; colIndex < iMatHyperDiffusion[dimIndex].cols(); colIndex++)
        {
          iMat_o[dimIndex](rowIndex,colIndex) = iMatHyperDiffusion[dimIndex](rowIndex,colIndex);
          lowerMat_o[dimIndex](rowIndex,colIndex) = lowerMatHyperDiffusion[dimIndex](rowIndex,colIndex);
          upperMat_o[dimIndex](rowIndex,colIndex) = upperMatHyperDiffusion[dimIndex](rowIndex,colIndex);
        }
      }
    }
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::getLowerReflectingBcMapping(unsigned dir,
        std::vector<unsigned>& nodeMap) const
  {
    throw Lucee::Except("CompletePolynomialElement::getLowerReflectingBcMapping: Not implemented");
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::getUpperReflectingBcMapping(unsigned dir,
        std::vector<unsigned>& nodeMap) const
  {
    throw Lucee::Except("CompletePolynomialElement::getUpperReflectingBcMapping: Not implemented");
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::getLowerFaceToInteriorMapping(unsigned dir,
        Lucee::Matrix<double>& faceToIntMap) const
  {
    throw Lucee::Except("CompletePolynomialElement::getLowerFaceToInteriorMapping: Not implemented");
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::extractFromField(const Lucee::Field<NDIM, double>& fld,
    std::vector<double>& data)
  {
    throw Lucee::Except("CompletePolynomialElement::extractFromField: Not implemented");
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::copyAllDataFromField(const Lucee::Field<NDIM, double>& fld,
    double *data)
  {
    throw Lucee::Except("CompletePolynomialElement::copyAllDataFromField: Not implemented!");
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::copyAllDataToField(const double *data, 
    Lucee::Field<NDIM, double>& fld)
  {
    throw Lucee::Except("CompletePolynomialElement::copyAllDataToField: Not implemented!");
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::evalBasis(double xc[NDIM], std::vector<double>& vals) const
  {
    // Create vector from xc
    Eigen::VectorXd nodeVec(NDIM);
    for (int componentIndex = 0; componentIndex < NDIM; componentIndex++)
      nodeVec(componentIndex) = xc[componentIndex];

    for (int basisIndex = 0; basisIndex < functionVector.size(); basisIndex++)
      vals[basisIndex] = evalPolynomial(functionVector[basisIndex], nodeVec);
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::setupMatrices()
  {
    // Get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->template getGrid<Lucee::StructuredGridBase<NDIM> >();
    // Get grid spacing
    for (int i = 0; i < NDIM; i++)
    {
      dq[i] = grid.getDx(i);
      dq2[i] = dq[i]*dq[i];
    }

    computeBasisFunctions(functionVector);

    // Compute gaussian quadrature weights and locations in 1-D
    numGaussPoints = (unsigned)((maxPower+1)/2.0 + 0.5);
    gaussPoints  = std::vector<double>(numGaussPoints);
    gaussWeights = std::vector<double>(numGaussPoints);
    legendre_set(numGaussPoints, &gaussPoints[0], &gaussWeights[0]);
    
    // Figure out how many gaussian points there are
    int totalVolumeGaussNodes = 1;

    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
      totalVolumeGaussNodes = gaussPoints.size()*totalVolumeGaussNodes;

    int totalSurfaceGaussNodes = totalVolumeGaussNodes/gaussPoints.size();
    
    gaussNodeList = Eigen::MatrixXd::Zero(totalVolumeGaussNodes, NDIM+1);
    functionEvaluations  = Eigen::MatrixXd(totalVolumeGaussNodes, this->getNumNodes());
    // 3D Array containing derivatives of basis functions (cols) evaluated at gaussian
    // integration locations (rows). Correspondance between column and gaussian node set
    // is kept track of in gaussNodeList
    blitz::Array<double,3> functionDEvaluations(totalVolumeGaussNodes, this->getNumNodes(),NDIM);
    
    gaussNodeListUpperSurf  = std::vector<Eigen::MatrixXd>(NDIM);
    gaussNodeListLowerSurf  = std::vector<Eigen::MatrixXd>(NDIM);
    upperSurfaceEvaluations = std::vector<Eigen::MatrixXd>(NDIM);
    lowerSurfaceEvaluations = std::vector<Eigen::MatrixXd>(NDIM);

    Eigen::MatrixXd gaussNodeListSurf = Eigen::MatrixXd::Zero(totalSurfaceGaussNodes, NDIM);

    if (NDIM > 1)
    {
      // TODO: roll surface and volume quadrature code into one function
      // since they do the same thing, but in diff dimensions
      // Compute surface gaussian quadrature locations
      int surfShape[NDIM-1];
      int surfIdx[NDIM-1];
      
      for (int dimIndex = 0; dimIndex < NDIM-1; dimIndex++)
        surfShape[dimIndex] = gaussPoints.size();

      Lucee::Region<NDIM-1, int> surfRegion(surfShape);

      Lucee::RowMajorSequencer<NDIM-1> surfSeq = RowMajorSequencer<NDIM-1>(surfRegion);
      Lucee::RowMajorIndexer<NDIM-1> surfIdxr = RowMajorIndexer<NDIM-1>(surfRegion);

      
      // Find all quadrature locations on a NDIM-1 surface
      while(surfSeq.step())
      {
        surfSeq.fillWithIndex(surfIdx);
        int nodeNumber = surfIdxr.getIndex(surfIdx);

        gaussNodeListSurf(nodeNumber, NDIM-1) = 1.0;
        for (int dimIndex = 0; dimIndex < NDIM-1; dimIndex++)
        {
          gaussNodeListSurf(nodeNumber, dimIndex) = gaussPoints[surfIdx[dimIndex]];
          gaussNodeListSurf(nodeNumber, NDIM-1)    *= gaussWeights[surfIdx[dimIndex]];
        }
      }
    }
    else
    {
      gaussNodeListSurf(0, 0) = 1.0;
    }

    // Evaluate quadrature points on all surfaces
    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      // Initialize matrices
      gaussNodeListUpperSurf[dimIndex]  = Eigen::MatrixXd::Zero(totalSurfaceGaussNodes, NDIM+1);
      gaussNodeListLowerSurf[dimIndex]  = Eigen::MatrixXd::Zero(totalSurfaceGaussNodes, NDIM+1);
      upperSurfaceEvaluations[dimIndex] = Eigen::MatrixXd::Zero(totalSurfaceGaussNodes, functionVector.size());
      lowerSurfaceEvaluations[dimIndex] = Eigen::MatrixXd::Zero(totalSurfaceGaussNodes, functionVector.size());
      
      Eigen::VectorXd gaussNodeVec = Eigen::VectorXd::Zero(NDIM);
      // Evaluate all basis functions at all upper and lower surfaces
      for (int nodeIndex = 0; nodeIndex < gaussNodeListSurf.rows(); nodeIndex++)
      {
        int rollingIndex = 0;
        
        for (int coordIndex = 0; coordIndex < NDIM; coordIndex++)
        {
          if (coordIndex == dimIndex)
            gaussNodeVec(dimIndex) = 1;
          else
          {
            gaussNodeVec(coordIndex) = gaussNodeListSurf(nodeIndex, rollingIndex);
            rollingIndex++;
          }
          gaussNodeListUpperSurf[dimIndex](nodeIndex, coordIndex) = gaussNodeVec(coordIndex);
        }
        gaussNodeListUpperSurf[dimIndex](nodeIndex, NDIM) = gaussNodeListSurf(nodeIndex, NDIM-1);

        // Evaluate all basis functions at this node location
        for (int basisIndex = 0; basisIndex < functionVector.size(); basisIndex++)
          upperSurfaceEvaluations[dimIndex](nodeIndex, basisIndex) = evalPolynomial(functionVector[basisIndex],gaussNodeVec);
        
        // Copy upper surface matrix into lower one first, then replace surface index
        gaussNodeListLowerSurf[dimIndex].row(nodeIndex) = gaussNodeListUpperSurf[dimIndex].row(nodeIndex);
        gaussNodeListLowerSurf[dimIndex](nodeIndex, dimIndex) = -1;
        gaussNodeVec(dimIndex) = -1;
        // Evaluate all basis functions at this node location
        for (int basisIndex = 0; basisIndex < functionVector.size(); basisIndex++)
          lowerSurfaceEvaluations[dimIndex](nodeIndex, basisIndex) = evalPolynomial(functionVector[basisIndex],gaussNodeVec);
      }
    }

    // Evaluate all basis functions and their derivatives in every direction at all gaussian volume nodes
    // First compute all volume gaussian quadrature locations
    int volShape[NDIM];
    int volIdx[NDIM];
    
    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
      volShape[dimIndex] = gaussPoints.size();

    Lucee::Region<NDIM, int> volRegion(volShape);

    // Needs col-major right now to help project 1-D data into 2-D element
    Lucee::ColMajorSequencer<NDIM> volSeq = ColMajorSequencer<NDIM>(volRegion);
    Lucee::ColMajorIndexer<NDIM> volIdxr = ColMajorIndexer<NDIM>(volRegion);

    // Find all quadrature locations on a NDIM volume
    while(volSeq.step())
    {
      volSeq.fillWithIndex(volIdx);
      int nodeNumber = volIdxr.getIndex(volIdx);

      gaussNodeList(nodeNumber, NDIM) = 1.0;
      Eigen::VectorXd gaussNodeVec = Eigen::VectorXd::Zero(NDIM);
      for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
      {
        gaussNodeVec(dimIndex)              = gaussPoints[volIdx[dimIndex]];
        gaussNodeList(nodeNumber, dimIndex) = gaussPoints[volIdx[dimIndex]];
        gaussNodeList(nodeNumber, NDIM)    *= gaussWeights[volIdx[dimIndex]];
      }

      for (int basisIndex = 0; basisIndex < functionVector.size(); basisIndex++)
      {
        functionEvaluations(nodeNumber, basisIndex) = evalPolynomial(functionVector[basisIndex],gaussNodeVec);
        
        for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
          functionDEvaluations(nodeNumber, basisIndex, dimIndex) = evalPolynomial(computePolynomialDerivative(functionVector[basisIndex],dimIndex),gaussNodeVec);
      }
    }

    // Resize+Initialize (most of) the output matrices we need
    resizeMatrices();

    // Call various functions to populate the matrices

    if (NDIM == 2)
    {
      // Comput moment matrices on reference element
      setupMomentMatrices();
      // Scale into physical space
      for(int pIndex = 0; pIndex < maxMoment+1; pIndex++)
        for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
          momMatrix[pIndex] *= 0.5*dq[dimIndex];
    }

    computeMass(refMass);

    computeStiffness(functionDEvaluations,refStiffness);

    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      computeFaceMass(dimIndex, refFaceMassLower[dimIndex], refFaceMassUpper[dimIndex]);
      computeGradStiffness(functionDEvaluations, dimIndex, refGradStiffness[dimIndex]);
    }

    // Scale the matrices computed on reference element into physical space
    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      refMass   *= 0.5*dq[dimIndex];
      refStiffness *= 0.5*dq[dimIndex];

      // Scale face-mass matrices
      for (int matrixIndex = 0; matrixIndex < NDIM; matrixIndex++)
      {
        if (dimIndex != matrixIndex)
        {
          refFaceMassUpper[matrixIndex] *= 0.5*dq[dimIndex];
          refFaceMassLower[matrixIndex] *= 0.5*dq[dimIndex];
          refGradStiffness[matrixIndex] *= 0.5*dq[dimIndex];
        }
      }
    }
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::resizeMatrices()
  {
    int numNodes = this->getNumNodes();

    refMass          = Eigen::MatrixXd::Zero(numNodes, numNodes);
    refStiffness     = Eigen::MatrixXd::Zero(numNodes, numNodes);

    refFaceMassLower = std::vector<Eigen::MatrixXd>(NDIM);
    refFaceMassUpper = std::vector<Eigen::MatrixXd>(NDIM);
    refGradStiffness = std::vector<Eigen::MatrixXd>(NDIM);

    for (unsigned dimIndex = 0; dimIndex < NDIM; dimIndex++)
    {
      refFaceMassLower[dimIndex] = Eigen::MatrixXd::Zero(numNodes, getNumSurfLowerNodes(dimIndex));
      refFaceMassUpper[dimIndex] = Eigen::MatrixXd::Zero(numNodes, getNumSurfUpperNodes(dimIndex));
      refGradStiffness[dimIndex] = Eigen::MatrixXd::Zero(numNodes, numNodes);
    }
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::setupBasisMatrix(Eigen::MatrixXi& basisMatrix)
  {
    int dataIndex = 0;
    int shape[NDIM];
    int idx[NDIM];

    for(int dimIndex = 0; dimIndex < NDIM; dimIndex++)
      shape[dimIndex] = polyOrder + 1;

    Lucee::Region<NDIM, int> polyRegion(shape);
    Lucee::RowMajorSequencer<NDIM> polySeq = RowMajorSequencer<NDIM>(polyRegion);

    while(polySeq.step())
    {
      polySeq.fillWithIndex(idx);

      // Compute total degree of monimial
      int totalDegree = 0;
      for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
        totalDegree += idx[dimIndex];
      
      if (totalDegree < polyOrder + 1)
      {
        for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
          basisMatrix(dataIndex, dimIndex) = idx[dimIndex];
        dataIndex++;
      }
    }
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::computeBasisFunctions(std::vector<blitz::Array<double,NDIM> >& functionVector)
  {
    // Matrix to represent basis monomials
    Eigen::MatrixXi basisList = Eigen::MatrixXi::Zero(this->getNumNodes(), NDIM);
    // Figure out basis monomials
    setupBasisMatrix(basisList);

    blitz::TinyVector<int, NDIM> polynomialShape(maxPower);
    blitz::TinyVector<int, NDIM> polynomialCoord;
    blitz::Array<double, NDIM> polynomialArray(polynomialShape);

    // Store each basis monomial into functionVector
    for (int basisIndex = 0; basisIndex < basisList.rows(); basisIndex++)
    {
      // Reset the polynomial3DArray
      polynomialArray = 0;
      // Figure out where this basis function will be stored
      polynomialCoord = 0;
      for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
        polynomialCoord(dimIndex) = basisList(basisIndex, dimIndex);

      polynomialArray(polynomialCoord) = 1.0;
      // Store the newly represented basis function into functionVector
      functionVector.push_back(polynomialArray.copy());
    }
  }

  template <unsigned NDIM>
  double
  CompletePolynomialElement<NDIM>::evalPolynomial(const blitz::Array<double,NDIM>& polyCoeffs, const VectorXd& nodeCoords) const
  {
    double totalSum = 0.0;
    double monomialTerm;
    int shape[NDIM];
    int idx[NDIM];
    blitz::TinyVector<int, NDIM> polyCoord;

    for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
      shape[dimIndex] = polyCoeffs.extent(dimIndex);

    Lucee::Region<NDIM, int> polyRegion(shape);
    Lucee::RowMajorSequencer<NDIM> polySeq = RowMajorSequencer<NDIM>(polyRegion);
    // Loop over each element of polyCoeffs
    while(polySeq.step())
    {
      polySeq.fillWithIndex(idx);
      for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
        polyCoord(dimIndex) = idx[dimIndex];
      
      if (polyCoeffs(polyCoord) != 0.0)
      {
        monomialTerm = polyCoeffs(polyCoord);
        // Compute monomial term associated with polyCoord
        for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
          monomialTerm *= pow(nodeCoords(dimIndex), idx[dimIndex]);
        totalSum += monomialTerm;
      }
    }

    return totalSum;
  }

  template <unsigned NDIM>
  blitz::Array<double, NDIM>
  CompletePolynomialElement<NDIM>::computePolynomialDerivative(const blitz::Array<double,NDIM>& poly, int dir)
  {
    int shape[NDIM];
    int idx[NDIM];
    blitz::TinyVector<int, NDIM> polyCoord;
    blitz::Array<double, NDIM> polyResult(poly.shape());
    polyResult = 0;

    for(int dimIndex = 0; dimIndex < NDIM; dimIndex++)
      shape[dimIndex] = poly.extent(dimIndex);

    Lucee::Region<NDIM, int> polyRegion(shape);

    Lucee::RowMajorSequencer<NDIM> polySeq = RowMajorSequencer<NDIM>(polyRegion);
    // Loop over each element of polyCoeffs
    while(polySeq.step())
    {
      polySeq.fillWithIndex(idx);
      for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
        polyCoord(dimIndex) = idx[dimIndex];
      
      if (polyCoord(dir) != 0 && poly(polyCoord) != 0.0)
      {
        double resultCoeff = polyCoord(dir)*poly(polyCoord);
        // Compute coordinate of derivative term
        polyCoord(dir) = polyCoord(dir) - 1;
        // Store new coefficient
        polyResult(polyCoord) += resultCoeff;
      }
    }

    return polyResult;
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::computeMass(Eigen::MatrixXd& resultMatrix)
  {
    for (int kIndex = 0; kIndex < resultMatrix.rows(); kIndex++)
    {
      for (int mIndex = 0; mIndex < resultMatrix.cols(); mIndex++)
      {
        // Reset integration result
        double integrationResult = 0.0;

        // Loop over 3d gauss points to evaluate integral using gaussian quadrature
        // result = weight(node) * f1(node) * f2(node)
        for (int gaussIndex = 0; gaussIndex < gaussNodeList.rows(); gaussIndex++)
          integrationResult += gaussNodeList(gaussIndex, NDIM)*functionEvaluations(gaussIndex, kIndex)*functionEvaluations(gaussIndex, mIndex);
        
        resultMatrix(kIndex, mIndex) = integrationResult;
      }
    }
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::computeFaceMass(int dir, Eigen::MatrixXd& lowerResultMatrix, Eigen::MatrixXd& upperResultMatrix)
  {
    std::vector<int> surfLowerNodeNums(lowerResultMatrix.cols(),0);
    std::vector<int> surfUpperNodeNums(upperResultMatrix.cols(),0);

    getSurfLowerNodeNums(dir,surfLowerNodeNums);
    getSurfUpperNodeNums(dir,surfUpperNodeNums);

    for (int kIndex = 0; kIndex < upperResultMatrix.rows(); kIndex++)
    {
      // Note that mIndex will be those of the ones in lower/upper node nums
      for (int mIndex = 0; mIndex < upperResultMatrix.cols(); mIndex++)
      {
        int lowerNodeNum = surfLowerNodeNums[mIndex];
        int upperNodeNum = surfUpperNodeNums[mIndex];
        double integrationResultL = 0.0;
        double integrationResultU = 0.0;
        
        for (int testIndex = 0; testIndex < gaussNodeListUpperSurf[dir].rows(); testIndex++)
        {
          integrationResultU += upperSurfaceEvaluations[dir](testIndex, upperNodeNum)*
            upperSurfaceEvaluations[dir](testIndex, kIndex)*gaussNodeListUpperSurf[dir](testIndex, NDIM);
          integrationResultL += lowerSurfaceEvaluations[dir](testIndex,lowerNodeNum)*
            lowerSurfaceEvaluations[dir](testIndex, kIndex)*gaussNodeListLowerSurf[dir](testIndex, NDIM);
        }

        upperResultMatrix(kIndex,mIndex) = integrationResultU;
        lowerResultMatrix(kIndex,mIndex) = integrationResultL;
      }
    }
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::computeStiffness(const blitz::Array<double, 3>& functionDerivative,
    Eigen::MatrixXd& resultMatrix)
  {
    for (int kIndex = 0; kIndex < resultMatrix.rows(); kIndex++)
    {
      for (int mIndex = 0; mIndex < resultMatrix.cols(); mIndex++)
      {
        // Reset integration result
        double integrationResult = 0.0;
        
        // Loop over 3d gauss points to evaluate integral using gaussian quadrature
        for (int gaussIndex = 0; gaussIndex < gaussNodeList.rows(); gaussIndex++)
          for (int dimIndex = 0; dimIndex < NDIM; dimIndex++)
            integrationResult += gaussNodeList(gaussIndex,NDIM)*functionDerivative(gaussIndex, kIndex, dimIndex)*functionDerivative(gaussIndex, mIndex, dimIndex)*4.0/dq2[dimIndex];
        
        resultMatrix(kIndex,mIndex) = integrationResult;
      }
    }
  }

  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::computeGradStiffness(const blitz::Array<double, 3>& functionDerivative,
    int dir, Eigen::MatrixXd& resultMatrix)
  {
    for (int kIndex = 0; kIndex < resultMatrix.rows(); kIndex++)
    {
      for (int mIndex = 0; mIndex < resultMatrix.cols(); mIndex++)
      {
        // Reset integration result
        double integrationResult = 0.0;

        // Loop over 3d gauss points to evaluate integral using gaussian quadrature
        for (int gaussIndex = 0; gaussIndex < gaussNodeList.rows(); gaussIndex++)
          integrationResult += gaussNodeList(gaussIndex,NDIM)*functionDerivative(gaussIndex, kIndex, dir)*functionEvaluations(gaussIndex, mIndex);

        resultMatrix(kIndex, mIndex) = integrationResult;
      }
    }
  }
 
  template <unsigned NDIM>
  void
  CompletePolynomialElement<NDIM>::setupMomentMatrices()
  {
    unsigned nodesPerSide = getNumSurfLowerNodes(1);
    std::vector<int> surfLowerNodeNums(nodesPerSide,1);
    getSurfLowerNodeNums(1,surfLowerNodeNums);

    maxMoment = 3;
    momMatrix = std::vector<Eigen::MatrixXd>(maxMoment+1);

    // Initialize Matrices
    for (int momentVal = 0; momentVal < maxMoment + 1; momentVal++)
      momMatrix[momentVal] = Eigen::MatrixXd::Zero(nodesPerSide, functionEvaluations.cols());

    // Evaluate Integral y^p*phi(x)psi(x,y)dA
    for (int j = 0; j < nodesPerSide; j++)
    {
      // 1-d basis function in this loop
      int lowerNodeNum = surfLowerNodeNums[j];
      for (int k = 0; k < functionEvaluations.cols(); k++)
      {
        for (int nodeIndex = 0; nodeIndex < gaussNodeList.rows(); nodeIndex++)
        {
          double gaussWeight  = gaussNodeList(nodeIndex,NDIM);
          double yCoord       = gaussNodeList(nodeIndex,1);
          unsigned rollingIndex = nodeIndex % numGaussPoints;

          double yPowerFactor = 1.0;
          for (int momentVal = 0; momentVal <= maxMoment; momentVal++)
          {
            if (momentVal != 0)
              yPowerFactor *= yCoord;
            momMatrix[momentVal](j,k) += gaussWeight*yPowerFactor*functionEvaluations(nodeIndex, k)*
              lowerSurfaceEvaluations[1](rollingIndex, lowerNodeNum);
          }
        }
      }
    }
  }

  template <unsigned NDIM>
  int
  CompletePolynomialElement<NDIM>::getCompletePolynomialDimension(int degree, int dimension) const
  {
    return factorial(dimension+degree)/(factorial(dimension)*factorial(degree));
  }

  template <unsigned NDIM>
  int
  CompletePolynomialElement<NDIM>::factorial(int n) const
  {
    return (n == 1 || n == 0) ? 1 : factorial(n-1)*n;
  }
  
  // instantiations
  template class CompletePolynomialElement<1>;
  template class CompletePolynomialElement<2>;
  template class CompletePolynomialElement<3>;
  template class CompletePolynomialElement<4>;
  template class CompletePolynomialElement<5>;
}
