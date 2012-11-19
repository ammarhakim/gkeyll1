/**
 * @file	LcSerendipityElement.cpp
 *
 * @brief       Interface class for a reference nodal finite element.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcExcept.h>
#include <LcGridIfc.h>
#include <LcSerendipityElement.h>

namespace Lucee
{
// set module name
  template <> const char *SerendipityElement<1>::id = "SerendipityElement";
  template <> const char *SerendipityElement<2>::id = "SerendipityElement";
  template <> const char *SerendipityElement<3>::id = "SerendipityElement";
  template <> const char *SerendipityElement<4>::id = "SerendipityElement";
  template <> const char *SerendipityElement<5>::id = "SerendipityElement";

  template <unsigned NDIM>
  SerendipityElement<NDIM>::SerendipityElement()
    : NodalFiniteElementIfc<NDIM>(1)
  {

  }
  
  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::NodalFiniteElementIfc<NDIM>::readInput(tbl);
    polyOrder = (unsigned) tbl.getNumber("polyOrder");
    // polyOrder is a misnomer for serendipity elements when we really want the degree
    basisDegree = polyOrder;

    if (basisDegree == 1)
      this->setNumNodes(8);
    else if (basisDegree == 2)
      this->setNumNodes(20);
    else if (basisDegree == 3)
      this->setNumNodes(32);
    else
    {
      Lucee::Except lce("SerendipityElement: Degree must be 1, 2, or 3.");
      lce << " Provided " << basisDegree << " instead";
      throw lce;
    }
  // Initialize matrices
    if (basisDegree == 1)
    {
      setupMatrices();
    }
    else if (basisDegree == 2)
    {
      setupMatrices();
    }
    else if (basisDegree == 3)
    {
      setupMatrices();
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getExclusiveNodeIndices(std::vector<int>& ndIds)
  {
    // Todo
    throw Lucee::Except("SerendipityElement::getExclusiveNodeIndices: Not implemented!");
  }

  template <unsigned NDIM>
  unsigned
  SerendipityElement<NDIM>::getNumSurfLowerNodes(unsigned dir) const
  {
    // Todo
    throw Lucee::Except("SerendipityElement::getNumSurfLowerNodes: Not implemented!");
    return 0;
  }

  template <unsigned NDIM>
  unsigned
  SerendipityElement<NDIM>::getNumSurfUpperNodes(unsigned dir) const
  {
    // Todo
    throw Lucee::Except("SerendipityElement::getNumSurfUpperNodes: Not implemented!");
    return 0;
  }

  template <unsigned NDIM>
  unsigned
  SerendipityElement<NDIM>::getNumGlobalNodes() const
  {
    throw Lucee::Except("SerendipityElement::getNumGlobalNodes: Not implemented!");
    return 0;
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getLocalToGlobal(std::vector<int>& lgMap) const
  {
    throw Lucee::Except("SerendipityElement::getLocalToGlobal: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getSurfLowerLocalToGlobal(unsigned dim,
    std::vector<int>& lgMap) const
  {
    throw Lucee::Except("SerendipityElement::getSurfLowerLocalToGlobal: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getSurfUpperLocalToGlobal(unsigned dim,
    std::vector<int>& lgMap) const
  {
    throw Lucee::Except("SerendipityElement::getSurfUpperLocalToGlobal: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getSurfLowerNodeNums(unsigned dir,
    std::vector<int>& nodeNum) const
  {
    // Todo
    throw Lucee::Except("SerendipityElement::getSurfLowerNodeNums: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getSurfUpperNodeNums(unsigned dir,
    std::vector<int>& nodeNum) const
  {
    // Todo
    throw Lucee::Except("SerendipityElement::getSurfUpperNodeNums: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getNodalCoordinates(Lucee::Matrix<double>& nodeCoords)
  {
    // Todo
    throw Lucee::Except("SerendipityElement::getNodalCoordinates: Not implemented!");    
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getWeights(std::vector<double>& w)
  {
    // Todo (2)
    throw Lucee::Except("SerendipityElement::getWeights: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getSurfUpperWeights(unsigned dir, std::vector<double>& w)
  {
    // Todo (2)
    throw Lucee::Except("SerendipityElement::getSurfUpperWeights: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getSurfLowerWeights(unsigned dir, std::vector<double>& w)
  {
    // Todo (2)
    throw Lucee::Except("SerendipityElement::getSurfLowerWeights: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getMassMatrix(Lucee::Matrix<double>& NjNk) const
  {
    // Todo
    throw Lucee::Except("SerendipityElement::getMassMatrix: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getLowerFaceMassMatrix(unsigned dir,
    Lucee::Matrix<double>& NjNk) const
  {
    // Todo
    throw Lucee::Except("SerendipityElement::getLowerFaceMassMatrix: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getUpperFaceMassMatrix(unsigned dir,
    Lucee::Matrix<double>& NjNk) const
  {
    // Todo
    throw Lucee::Except("SerendipityElement::getUpperFaceMassMatrix: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getStiffnessMatrix(Lucee::Matrix<double>& DNjDNk) const
  {
    throw Lucee::Except("SerendipityElement::getStiffnessMatrix: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getGradStiffnessMatrix(
    unsigned dir, Lucee::Matrix<double>& DNjNk) const
  {
    // Todo
    throw Lucee::Except("SerendipityElement::getGradStiffnessMatrix: Not implemented!");
  }

  template <unsigned NDIM>
  unsigned
  SerendipityElement<NDIM>::getNumGaussNodes() const
  {
    throw Lucee::Except("SerendipityElement::getNumGaussNodes: Not implemented!");
  }

  template <unsigned NDIM>
  unsigned
  SerendipityElement<NDIM>::getNumSurfGaussNodes() const
  {
    throw Lucee::Except("SerendipityElement::getNumSurfGaussNodes: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getGaussQuadData(Lucee::Matrix<double>& interpMat,
    Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const
  {
    throw Lucee::Except("SerendipityElement::getGaussQuadData: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getSurfLowerGaussQuadData(unsigned dir, Lucee::Matrix<double>& interpMat,
    Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const
  {
    throw Lucee::Except("SerendipityElement::getSurfLowerGaussQuadData: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getSurfUpperGaussQuadData(unsigned dir, Lucee::Matrix<double>& interpMat,
    Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const
  {
    throw Lucee::Except("SerendipityElement::getSurfUpperGaussQuadData: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getMomentMatrix(unsigned p, Lucee::Matrix<double>& momMat) const
  {
    throw Lucee::Except("SerendipityElement::getMomentMatrix: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::extractFromField(const Lucee::Field<NDIM, double>& fld,
    std::vector<double>& data)
  {
    throw Lucee::Except("SerendipityElement::extractFromField: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::copyAllDataFromField(const Lucee::Field<NDIM, double>& fld,
    double *data)
  {
    throw Lucee::Except("SerendipityElement::copyAllDataFromField: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::copyAllDataToField(const double *data, 
    Lucee::Field<NDIM, double>& fld)
  {
    throw Lucee::Except("SerendipityElement::copyAllDataToField: Not implemented!");
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::setupMatrices()
  {
    // Get hold of grid
    const Lucee::StructuredGridBase<2>& grid = this->template getGrid<Lucee::StructuredGridBase<2> >();
    // Get grid spacing (this is assumed to be uniform for now)
    double dx = grid.getDx(0);
    double dy = grid.getDx(1);
    // TODO: change when we test in 3d
    double dz = grid.getDx(1);
    double dx2 = dx*dx;
    double dy2 = dy*dy;
    double dz2 = dz*dz;
    // Populate nodeList according to basisDegree
    // TODO: Give Matrix name in argument for readability
    nodeList = MatrixXd(this->getNumNodes(),3);
    setupNodeListMatrix();
    // Populate basisList according to basisDegree
    basisList = MatrixXi(this->getNumNodes(),3);
    setupBasisMatrix();
    // Compute coefficients for basis functions
    MatrixXd coeffMatrix(this->getNumNodes(),this->getNumNodes());
    double coeffProduct;

    // TODO: Write a function that only sets the values of functionVector

    // Compute monomial terms evaluated at each nodal coordinate
    for (unsigned nodeIndex = 0; nodeIndex < nodeList.rows(); nodeIndex++)
    {
      for (unsigned basisIndex = 0; basisIndex < basisList.rows(); basisIndex++)
      {
        coeffProduct = 1.0;
        for (unsigned dimIndex = 0; dimIndex < basisList.cols(); dimIndex++)
        {
          coeffProduct = coeffProduct*pow(nodeList(nodeIndex,dimIndex),basisList(basisIndex,dimIndex));
        }
        coeffMatrix(nodeIndex,basisIndex) = coeffProduct;
      }
    }

    // Each column of invCoeffMatrix will be coefficient of basis monomial in basisList
    MatrixXd invCoeffMatrix = coeffMatrix.inverse();

    // Check that basis functions evaluate to 1 at each node and 0 otherwise:
    MatrixXd evalMatrix(this->getNumNodes(), this->getNumNodes());
    // Compute maximum polynomial power we want (no transformations)
    maxPower = 2*(polyOrder + 3 - 1); 
    blitz::Array<double,3> polynomial3DArray(maxPower,maxPower,maxPower);
    std::vector<blitz::Array<double,3> > functionVector;

    int xPow,yPow,zPow;

    // Loop over the coefficient matrix (each column represents a basis function)
    for (unsigned functionIndex = 0; functionIndex < invCoeffMatrix.cols(); functionIndex++)
    {
      // Reset the polynomial3DArray
      polynomial3DArray = 0;
      // Loop over the monomial terms that make up the basis function
      for (unsigned monomialIndex = 0; monomialIndex < invCoeffMatrix.rows(); monomialIndex++)
      {
        // Store the basis function in polynomial3DArray by storing its coefficient at
        // a location specified by the monimial term powers
        xPow = basisList(monomialIndex,0);
        yPow = basisList(monomialIndex,1);
        zPow = basisList(monomialIndex,2);
        polynomial3DArray(xPow,yPow,zPow) = polynomial3DArray(xPow,yPow,zPow) + invCoeffMatrix(monomialIndex,functionIndex);
      }
      // Store the newly represented basis function into functionVector
      functionVector.push_back(polynomial3DArray.copy());
    }
    
    // Below is for testing purposes only!
    /**for (unsigned nodeIndex = 0; nodeIndex < nodeList.rows(); nodeIndex++)
    {
      // Loop over the basis functions
      for (unsigned basisIndex = 0; basisIndex < invCoeffMatrix.cols(); basisIndex++)
      {
          // Evaluate each basis function at the node we are at to check for correct behavior
          evalMatrix(nodeIndex,basisIndex) = evalPolynomial(functionVector[basisIndex],nodeList.row(nodeIndex));
      }
    }*/

    // Compute gaussian quadrature weights and locations in 1-D
    // TODO: should these be at the class scope?
    unsigned numGaussPoints = (unsigned)((maxPower-1)/2.0 + 0.5);
    gaussPoints  = std::vector<double>(numGaussPoints);
    gaussWeights = std::vector<double>(numGaussPoints);
    legendre_set(numGaussPoints, &gaussPoints[0], &gaussWeights[0]);
    
    // Resize the various matrices now we know the node count
    // TODO: check if this notation is necessary
    refNjNk.resize(this->getNumNodes(),this->getNumNodes());
    refDNjNk_0.resize(this->getNumNodes(),this->getNumNodes());
    refDNjNk_1.resize(this->getNumNodes(),this->getNumNodes());
    refDNjNk_2.resize(this->getNumNodes(),this->getNumNodes());
    refFaceNjNk_xl.resize(this->getNumNodes(),this->getNumNodes());
    refFaceNjNk_xu.resize(this->getNumNodes(),this->getNumNodes());
    refFaceNjNk_yl.resize(this->getNumNodes(),this->getNumNodes());
    refFaceNjNk_yu.resize(this->getNumNodes(),this->getNumNodes());
    refFaceNjNk_zl.resize(this->getNumNodes(),this->getNumNodes());
    refFaceNjNk_zu.resize(this->getNumNodes(),this->getNumNodes());
    // Call various functions to populate the matrices
    refNjNk    = computeMass(functionVector);
    computeFaceMass(functionVector,0,refFaceNjNk_xl,refFaceNjNk_xu);
    computeFaceMass(functionVector,1,refFaceNjNk_yl,refFaceNjNk_yu);
    computeFaceMass(functionVector,2,refFaceNjNk_zl,refFaceNjNk_zu);
    refDNjDNk  = computeStiffness(functionVector);
    refDNjNk_0 = computeGradStiffness(functionVector,0);
    refDNjNk_1 = computeGradStiffness(functionVector,1);
    refDNjNk_2 = computeGradStiffness(functionVector,2);
    // Scale the matrices computed on reference element into physical space
    // TODO: verify correctness of this
    refNjNk        *= 0.5*dx*0.5*dy*0.5*dz;
    refFaceNjNk_xl *= 0.5*dy*0.5*dz;
    refFaceNjNk_xu *= 0.5*dy*0.5*dz;
    refFaceNjNk_yl *= 0.5*dx*0.5*dz;
    refFaceNjNk_yu *= 0.5*dx*0.5*dz;
    refFaceNjNk_zl *= 0.5*dx*0.5*dy;
    refFaceNjNk_zu *= 0.5*dx*0.5*dy;
    refDNjDNk      *= 0.5*dx*0.5*dy*0.5*dz;
    refDNjNk_0     *= 0.5*dx*0.5*dy*0.5*dz;
    refDNjNk_1     *= 0.5*dx*0.5*dy*0.5*dz;
    refDNjNk_2     *= 0.5*dx*0.5*dy*0.5*dz;
  }
 
  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::setupNodeListMatrix()
  {
    if (basisDegree == 1) {
      nodeList << -1,-1,-1,
                   1,-1,-1,
                   1,1,-1,
                   -1,1,-1,
                   -1,-1,1,
                   1,-1,1,
                   1,1,1,
                   -1,1,1;
    }
    else if (basisDegree == 2)
    {
      nodeList << -1,-1,-1,
                  0,-1,-1,
                  1,-1,-1,
                  1,0,-1,
                  1,1,-1,
                  0,1,-1,
                  -1,1,-1,
                  -1,0,-1,// End of Bottom Layer
                  -1,-1,0,
                  1,-1,0,
                  1,1,0,
                  -1,1,0, // End of Middle Layer
                  -1,-1,1,
                  0,-1,1,
                  1,-1,1,
                  1,0,1,
                  1,1,1,
                  0,1,1,
                  -1,1,1,
                  -1,0,1; // End of Top Layer
    }
    else if (basisDegree == 3)
    {
      nodeList << -1,-1,-1,
                  -1/3.0,-1,-1,
                  1/3.0,-1,-1,
                  1,-1,-1,
                  1,-1/3.0,-1,
                  1,1/3.0,-1,
                  1,1,-1,
                  1/3.0,1,-1,
                  -1/3.0,1,-1,
                  -1,1,-1,
                  -1,1/3.0,-1,
                  -1,-1/3.0,-1,// End of Bottom Layer
                  -1,-1,-1/3.0,
                  1,-1,-1/3.0,
                  1,1,-1/3.0,
                  -1,1,-1/3.0, // End of 2nd Layer
                  -1,-1,1/3.0,
                  1,-1,1/3.0,
                  1,1,1/3.0,
                  -1,1,1/3.0,  // End of 3rd Layer
                  -1,-1,1,
                  -1/3.0,-1,1,
                  1/3.0,-1,1,
                  1,-1,1,
                  1,-1/3.0,1,
                  1,1/3.0,1,
                  1,1,1,
                  1/3.0,1,1,
                  -1/3.0,1,1,
                  -1,1,1,
                  -1,1/3.0,1,
                  -1,-1/3.0,1; // End of Top Layer
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::setupBasisMatrix()
  {
    // Superlinear degree
    unsigned superDegree;
    unsigned dataIndex = 0;
    for (unsigned zIndex = 0; zIndex <= basisDegree; zIndex++)
    {
      for (unsigned yIndex = 0; yIndex <= basisDegree; yIndex++)
      {
        for (unsigned xIndex = 0; xIndex <= basisDegree; xIndex++)
        {
          superDegree = (xIndex>1)*xIndex + (yIndex>1)*yIndex + (zIndex>1)*zIndex;
          
          if (superDegree < basisDegree + 1)
          {
            basisList.row(dataIndex) = Vector3i(xIndex,yIndex,zIndex);
            dataIndex = dataIndex + 1;
          }
        }
      }
    }
  }

  template <unsigned NDIM>
  double
  SerendipityElement<NDIM>::evalPolynomial(blitz::Array<double,3> polyCoeffs, VectorXd nodeCoords)
  {
    double totalSum = 0.0;
    double monomialTerm;

    for (int xIndex = 0; xIndex < polyCoeffs.extent(0); xIndex++)
    {
      for (int yIndex = 0; yIndex < polyCoeffs.extent(1); yIndex++)
      {
        for (int zIndex = 0; zIndex < polyCoeffs.extent(2); zIndex++)
        {
          // Maybe add a check to see if polyCoeffs is 0 so we don't have to compute pow's
          monomialTerm = pow(nodeCoords(0),xIndex)*pow(nodeCoords(1),yIndex)*pow(nodeCoords(2),zIndex);
          totalSum += polyCoeffs(xIndex,yIndex,zIndex)*monomialTerm;
        }
      }
    }

    return totalSum;
  }
  
  template <unsigned NDIM>
  blitz::Array<double,3>
  SerendipityElement<NDIM>::computePolynomialProduct(blitz::Array<double,3> poly1, blitz::Array<double,3> poly2)
  {
    int x3Index;
    int y3Index;
    int z3Index;

    blitz::Array<double,3> poly3(poly1.extent(0),poly1.extent(1),poly1.extent(2));
    poly3 = 0;
    
    // Loop over each element of poly1, multiplying it by every element of poly2...
    for (int x1Index = 0; x1Index < poly1.extent(0); x1Index++)
    {
      for (int y1Index = 0; y1Index < poly1.extent(1); y1Index++)
      {
        for (int z1Index = 0; z1Index < poly1.extent(2); z1Index++)
        {
          for (int x2Index = 0; x2Index < poly2.extent(0); x2Index++)
          {
            for (int y2Index = 0; y2Index < poly2.extent(1); y2Index++)
            {
              for (int z2Index = 0; z2Index < poly2.extent(2); z2Index++)
              {
                x3Index = x1Index + x2Index;
                y3Index = y1Index + y2Index;
                z3Index = z1Index + z2Index;
                
                if ((x3Index > poly3.extent(0)) || (y3Index > poly3.extent(1)) || (z3Index > poly3.extent(2)))
                  continue;

                poly3(x3Index,y3Index,z3Index) += poly1(x1Index,y1Index,z1Index)*poly2(x2Index,y2Index,z2Index);
              }
            }
          }
        }
      }
    }

    return poly3;
  }

  template <unsigned NDIM>
  blitz::Array<double,3>
  SerendipityElement<NDIM>::computePolynomialDerivative(blitz::Array<double,3> poly, unsigned dir)
  {
    int xResultIndex;
    int yResultIndex;
    int zResultIndex;
    double resultCoefficient;

    blitz::Array<double,3> polyResult(poly.extent(0),poly.extent(1),poly.extent(2));
    polyResult = 0;
    
    // Loop over each element of poly1, multiplying it by every element of poly2...
    for (int xIndex = 0; xIndex < poly.extent(0); xIndex++)
    {
      for (int yIndex = 0; yIndex < poly.extent(1); yIndex++)
      {
        for (int zIndex = 0; zIndex < poly.extent(2); zIndex++)
        {
          if (dir == 0)
          {
            // Compute new location of coefficient
            if (xIndex == 0)
              continue;
            else {
              xResultIndex = xIndex - 1;
              yResultIndex = yIndex;
              zResultIndex = zIndex;
            }
            // Compute and store new coefficient
            polyResult(xResultIndex,yResultIndex,zResultIndex) += xIndex*poly(xIndex,yIndex,zIndex);
          }
          else if (dir == 1)
          {
            // Compute new location of coefficient
            if (yIndex == 0)
              continue;
            else {
              xResultIndex = xIndex;
              yResultIndex = yIndex - 1;
              zResultIndex = zIndex;
            }
            // Compute and store new coefficient
            polyResult(xResultIndex,yResultIndex,zResultIndex) += yIndex*poly(xIndex,yIndex,zIndex);
          }
          else if (dir == 2)
          {
            // Compute new location of coefficient
            if (zIndex == 0)
              continue;
            else {
              xResultIndex = xIndex;
              yResultIndex = yIndex;
              zResultIndex = zIndex-1;
            }
            // Compute and store new coefficient
            polyResult(xResultIndex,yResultIndex,zResultIndex) += zIndex*poly(xIndex,yIndex,zIndex);
          }
        }
      }
    }

    return polyResult;
  }

  template <unsigned NDIM>
  blitz::Array<double,2>
  SerendipityElement<NDIM>::computeMass(std::vector<blitz::Array<double,3> > functionVector)
  {
    blitz::Array<double,3> polyProduct(maxPower,maxPower,maxPower);
    blitz::Array<double,2> resultArray(this->getNumNodes(),this->getNumNodes());
    VectorXd gaussNodeVec(3);

    resultArray = 0;
    double integrationResult;
    double totalWeight;

    for (unsigned kIndex = 0; kIndex < resultArray.rows(); kIndex++)
    {
      for (unsigned mIndex = 0; mIndex < resultArray.cols(); mIndex++)
      {
        // Reset polyProduct
        polyProduct = 0;
        // Calculate polynomial product of two basis functions
        polyProduct = computePolynomialProduct(functionVector[kIndex],functionVector[mIndex]);

        // Reset integration result
        integrationResult = 0.0;
        // Integrate using gaussian quadrature
        for (unsigned xNodeIndex = 0; xNodeIndex < gaussPoints.size(); xNodeIndex++)
        {
          for (unsigned yNodeIndex = 0; yNodeIndex < gaussPoints.size(); yNodeIndex++)
          {
            for (unsigned zNodeIndex = 0; zNodeIndex < gaussPoints.size(); zNodeIndex++)
            {
              gaussNodeVec(0)    = gaussPoints[xNodeIndex];
              gaussNodeVec(1)    = gaussPoints[yNodeIndex];
              gaussNodeVec(2)    = gaussPoints[zNodeIndex];
              totalWeight        = gaussWeights[xNodeIndex]*gaussWeights[yNodeIndex]*gaussWeights[zNodeIndex];
              integrationResult += totalWeight*evalPolynomial(polyProduct,gaussNodeVec);
            }
          }
        }

        resultArray((int) kIndex,(int) mIndex) = integrationResult;
      }
    }

    return resultArray;
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::computeFaceMass(std::vector<blitz::Array<double,3> > functionVector, unsigned dir,
                                            blitz::Array<double,2> &lowerResultArray, 
                                            blitz::Array<double,2> &upperResultArray)
  {
    blitz::Array<double,3> polyProduct(maxPower,maxPower,maxPower);
    VectorXd gaussNodeVec(3);

    upperResultArray = 0;
    lowerResultArray = 0;
    double integrationResultU;
    double integrationResultL;
    double totalWeight;

    for (unsigned kIndex = 0; kIndex < upperResultArray.rows(); kIndex++)
    {
      for (unsigned mIndex = 0; mIndex < upperResultArray.cols(); mIndex++)
      {
        // Reset polyProduct
        polyProduct = 0;
        // Calculate polynomial product of two basis functions
        polyProduct = computePolynomialProduct(functionVector[kIndex],functionVector[mIndex]);

        // Reset integration result
        integrationResultL = 0.0;
        integrationResultU = 0.0;
        // Integrate using gaussian quadrature
        for (unsigned dir1NodeIndex = 0; dir1NodeIndex < gaussPoints.size(); dir1NodeIndex++)
        {
          for (unsigned dir2NodeIndex = 0; dir2NodeIndex < gaussPoints.size(); dir2NodeIndex++)
          {
            if (dir == 0)
            {
              gaussNodeVec(0) = 1;
              gaussNodeVec(1) = gaussPoints[dir1NodeIndex];
              gaussNodeVec(2) = gaussPoints[dir2NodeIndex];
              totalWeight = gaussWeights[dir1NodeIndex]*gaussWeights[dir2NodeIndex];
            }
            else if (dir == 1)
            {
              gaussNodeVec(0) = gaussPoints[dir1NodeIndex];
              gaussNodeVec(1) = 1;
              gaussNodeVec(2) = gaussPoints[dir2NodeIndex];
              totalWeight = gaussWeights[dir1NodeIndex]*gaussWeights[dir2NodeIndex];
            }
            else if (dir == 2)
            {
              gaussNodeVec(0) = gaussPoints[dir1NodeIndex];
              gaussNodeVec(1) = gaussPoints[dir2NodeIndex];
              gaussNodeVec(2) = 1;
              totalWeight = gaussWeights[dir1NodeIndex]*gaussWeights[dir2NodeIndex];
            }
            integrationResultU += totalWeight*evalPolynomial(polyProduct,gaussNodeVec);
            // Replace the appropriate coordinate to evaluation location
            gaussNodeVec(dir) = -1;
            integrationResultL += totalWeight*evalPolynomial(polyProduct,gaussNodeVec);
          }
        }

        upperResultArray((int) kIndex,(int) mIndex) = integrationResultU;
        lowerResultArray((int) kIndex,(int) mIndex) = integrationResultL;
      }
    }
  }

  template <unsigned NDIM>
  blitz::Array<double,2>
  SerendipityElement<NDIM>::computeStiffness(std::vector<blitz::Array<double,3> > functionVector)
  {
    blitz::Array<double,2> resultArray(this->getNumNodes(),this->getNumNodes());
    blitz::Array<double,3> polyProduct(maxPower,maxPower,maxPower);
    blitz::Array<double,3> functionDx2(maxPower,maxPower,maxPower);
    blitz::Array<double,3> functionDy2(maxPower,maxPower,maxPower);
    blitz::Array<double,3> functionDz2(maxPower,maxPower,maxPower);
    VectorXd gaussNodeVec(3);

    resultArray = 0;
    double integrationResult;
    double totalWeight;
    // Compute stiffness matrix (K)
    
    
    for (unsigned kIndex = 0; kIndex < resultArray.rows(); kIndex++)
    {
      for (unsigned mIndex = 0; mIndex < resultArray.cols(); mIndex++)
      {
        // Reset polyProduct
        polyProduct = 0;
        functionDx2 = 0;
        functionDy2 = 0;
        functionDz2 = 0;
        // Calculate derivatives of basis functions
        functionDx2 = computePolynomialProduct(computePolynomialDerivative(functionVector[kIndex],0),
                                               computePolynomialDerivative(functionVector[mIndex],0));
        functionDy2 = computePolynomialProduct(computePolynomialDerivative(functionVector[kIndex],1),
                                               computePolynomialDerivative(functionVector[mIndex],1));
        functionDz2 = computePolynomialProduct(computePolynomialDerivative(functionVector[kIndex],2),
                                               computePolynomialDerivative(functionVector[mIndex],2));
        // Calculate polynomial product of two basis functions
        polyProduct = functionDx2 + functionDy2 + functionDz2;

        // Reset integration result
        integrationResult = 0.0;
        // Integrate using gaussian quadrature
        for (unsigned xNodeIndex = 0; xNodeIndex < gaussPoints.size(); xNodeIndex++)
        {
          for (unsigned yNodeIndex = 0; yNodeIndex < gaussPoints.size(); yNodeIndex++)
          {
            for (unsigned zNodeIndex = 0; zNodeIndex < gaussPoints.size(); zNodeIndex++)
            {
            gaussNodeVec(0)    = gaussPoints[xNodeIndex];
            gaussNodeVec(1)    = gaussPoints[yNodeIndex];
            gaussNodeVec(2)    = gaussPoints[zNodeIndex];
            totalWeight        = gaussWeights[xNodeIndex]*gaussWeights[yNodeIndex]*gaussWeights[zNodeIndex];
            integrationResult += totalWeight*evalPolynomial(polyProduct,gaussNodeVec);
            }
          }
        }

        resultArray((int) kIndex,(int) mIndex) = integrationResult;
      }
    }
    
    return resultArray;
  }

  template <unsigned NDIM>
  blitz::Array<double,2>
  SerendipityElement<NDIM>::computeGradStiffness(std::vector<blitz::Array<double,3> > functionVector, 
                                                 unsigned dir)
  {
    blitz::Array<double,3> polyProduct(maxPower,maxPower,maxPower);
    blitz::Array<double,2> resultArray(this->getNumNodes(),this->getNumNodes());
    blitz::Array<double,3> functionD(maxPower,maxPower,maxPower);
    VectorXd gaussNodeVec(3);

    resultArray = 0;
    double integrationResult;
    double totalWeight;

    for (unsigned kIndex = 0; kIndex < resultArray.rows(); kIndex++)
    {
      for (unsigned mIndex = 0; mIndex < resultArray.cols(); mIndex++)
      {
        // Reset polyProduct
        polyProduct = 0;
        functionD = 0;
        // Calculate derivatives of basis functions
        functionD = computePolynomialDerivative(functionVector[kIndex],dir);
        // Calculate polynomial product of two basis functions
        polyProduct = computePolynomialProduct(functionD,functionVector[mIndex]);

        // Reset integration result
        integrationResult = 0.0;
        // Integrate using gaussian quadrature
        for (unsigned xNodeIndex = 0; xNodeIndex < gaussPoints.size(); xNodeIndex++)
        {
          for (unsigned yNodeIndex = 0; yNodeIndex < gaussPoints.size(); yNodeIndex++)
          {
            for (unsigned zNodeIndex = 0; zNodeIndex < gaussPoints.size(); zNodeIndex++)
            {
              gaussNodeVec(0)    = gaussPoints[xNodeIndex];
              gaussNodeVec(1)    = gaussPoints[yNodeIndex];
              gaussNodeVec(2)    = gaussPoints[zNodeIndex];
              totalWeight        = gaussWeights[xNodeIndex]*gaussWeights[yNodeIndex]*gaussWeights[zNodeIndex];
              integrationResult += totalWeight*evalPolynomial(polyProduct,gaussNodeVec);
            }
          }
        }

        resultArray((int) kIndex,(int) mIndex) = integrationResult;
      }
    }

    return resultArray;
  }
  
// instantiations
  template class SerendipityElement<1>;
  template class SerendipityElement<2>;
  template class SerendipityElement<3>;
  template class SerendipityElement<4>;
  template class SerendipityElement<5>;
}
