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
  using namespace Eigen;

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

    if (NDIM == 3)
    {
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

      setupMatrices();
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getExclusiveNodeIndices(std::vector<int>& ndIds)
  {
    throw Lucee::Except("SerendipityElement::getExclusiveNodeIndices: Not implemented!");
  }

  template <unsigned NDIM>
  unsigned
  SerendipityElement<NDIM>::getNumSurfLowerNodes(unsigned dir) const
  {
    if (basisDegree == 1)
    {
      return 4;
    }
    else if (basisDegree == 2)
    {
      return 8;
    }
    else if (basisDegree == 3)
    {
      return 12;
    }
  }

  template <unsigned NDIM>
  unsigned
  SerendipityElement<NDIM>::getNumSurfUpperNodes(unsigned dir) const
  {
    if (basisDegree == 1)
    {
      return 4;
    }
    else if (basisDegree == 2)
    {
      return 8;
    }
    else if (basisDegree == 3)
    {
      return 12;
    }
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
    if (basisDegree == 1)
    {
      if (dir == 0)
      {
        nodeNum[0] = 0;
        nodeNum[1] = 3;
        nodeNum[2] = 4;
        nodeNum[3] = 7;
      }
      else if (dir == 1)
      {
        nodeNum[0] = 0;
        nodeNum[1] = 1;
        nodeNum[2] = 4;
        nodeNum[3] = 5;
      }
      else if (dir == 2)
      {
        nodeNum[0] = 0;
        nodeNum[1] = 1;
        nodeNum[2] = 3;
        nodeNum[3] = 2;
      }
    }
    else if (basisDegree == 2)
    {
      if (dir == 0)
      {
        nodeNum[0] = 0;
        nodeNum[1] = 7;
        nodeNum[2] = 6;
        nodeNum[3] = 8;
        nodeNum[4] = 11;
        nodeNum[5] = 12;
        nodeNum[6] = 19;
        nodeNum[7] = 18;
      }
      else if (dir == 1)
      {
        nodeNum[0] = 0;
        nodeNum[1] = 1;
        nodeNum[2] = 2;
        nodeNum[3] = 8;
        nodeNum[4] = 9;
        nodeNum[5] = 12;
        nodeNum[6] = 13;
        nodeNum[7] = 14;
      }
      else if (dir == 2)
      {
        nodeNum[0] = 0;
        nodeNum[1] = 1;
        nodeNum[2] = 2;
        nodeNum[3] = 7;
        nodeNum[4] = 3;
        nodeNum[5] = 6;
        nodeNum[6] = 5;
        nodeNum[7] = 4;
      }
    }
    else if (basisDegree == 3)
    {
      if (dir == 0)
      {
        nodeNum[0] = 0;
        nodeNum[1] = 11;
        nodeNum[2] = 10;
        nodeNum[3] = 9;
        nodeNum[4] = 12;
        nodeNum[5] = 15;
        nodeNum[6] = 16;
        nodeNum[7] = 19;
        nodeNum[8] = 20;
        nodeNum[9] = 31;
        nodeNum[10] = 30;
        nodeNum[11] = 29;
      }
      else if (dir == 1)
      {
        nodeNum[0] = 0;
        nodeNum[1] = 1;
        nodeNum[2] = 2;
        nodeNum[3] = 3;
        nodeNum[4] = 12;
        nodeNum[5] = 13;
        nodeNum[6] = 16;
        nodeNum[7] = 17;
        nodeNum[8] = 20;
        nodeNum[9] = 21;
        nodeNum[10] = 22;
        nodeNum[11] = 23;
      }
      else if (dir == 2)
      {
        nodeNum[0] = 0;
        nodeNum[1] = 1;
        nodeNum[2] = 2;
        nodeNum[3] = 3;
        nodeNum[4] = 11;
        nodeNum[5] = 4;
        nodeNum[6] = 10;
        nodeNum[7] = 5;
        nodeNum[8] = 9;
        nodeNum[9] = 8;
        nodeNum[10] = 7;
        nodeNum[11] = 6;
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getSurfUpperNodeNums(unsigned dir,
    std::vector<int>& nodeNum) const
  {
    if (basisDegree == 1)
    {
      if (dir == 0)
      {
        nodeNum[0] = 1;
        nodeNum[1] = 2;
        nodeNum[2] = 5;
        nodeNum[3] = 6;
      }
      else if (dir == 1)
      {
        nodeNum[0] = 3;
        nodeNum[1] = 2;
        nodeNum[2] = 7;
        nodeNum[3] = 6;
      }
      else if (dir == 2)
      {
        nodeNum[0] = 4;
        nodeNum[1] = 5;
        nodeNum[2] = 7;
        nodeNum[3] = 6;
      }
    }
    else if (basisDegree == 2)
    {
      if (dir == 0)
      {
        nodeNum[0] = 2;
        nodeNum[1] = 3;
        nodeNum[2] = 4;
        nodeNum[3] = 9;
        nodeNum[4] = 10;
        nodeNum[5] = 14;
        nodeNum[6] = 15;
        nodeNum[7] = 16;
      }
      else if (dir == 1)
      {
        nodeNum[0] = 6;
        nodeNum[1] = 5;
        nodeNum[2] = 4;
        nodeNum[3] = 11;
        nodeNum[4] = 10;
        nodeNum[5] = 18;
        nodeNum[6] = 17;
        nodeNum[7] = 16;
      }
      else if (dir == 2)
      {
        nodeNum[0] = 12;
        nodeNum[1] = 13;
        nodeNum[2] = 14;
        nodeNum[3] = 19;
        nodeNum[4] = 15;
        nodeNum[5] = 18;
        nodeNum[6] = 17;
        nodeNum[7] = 16;
      }
    }
    else if (basisDegree == 3)
    {
      if (dir == 0)
      {
        nodeNum[0] = 3;
        nodeNum[1] = 4;
        nodeNum[2] = 5;
        nodeNum[3] = 6;
        nodeNum[4] = 13;
        nodeNum[5] = 14;
        nodeNum[6] = 17;
        nodeNum[7] = 18;
        nodeNum[8] = 23;
        nodeNum[9] = 24;
        nodeNum[10] = 25;
        nodeNum[11] = 26;
      }
      else if (dir == 1)
      {
        nodeNum[0] = 9;
        nodeNum[1] = 8;
        nodeNum[2] = 7;
        nodeNum[3] = 6;
        nodeNum[4] = 15;
        nodeNum[5] = 14;
        nodeNum[6] = 19;
        nodeNum[7] = 18;
        nodeNum[8] = 29;
        nodeNum[9] = 28;
        nodeNum[10] = 27;
        nodeNum[11] = 26;
      }
      else if (dir == 2)
      {
        nodeNum[0] = 20;
        nodeNum[1] = 21;
        nodeNum[2] = 22;
        nodeNum[3] = 23;
        nodeNum[4] = 31;
        nodeNum[5] = 24;
        nodeNum[6] = 30;
        nodeNum[7] = 25;
        nodeNum[8] = 29;
        nodeNum[9] = 28;
        nodeNum[10] = 27;
        nodeNum[11] = 26;
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getNodalCoordinates(Lucee::Matrix<double>& nodeCoords)
  {
    const Lucee::StructuredGridBase<3>& grid 
      = this->template getGrid<Lucee::StructuredGridBase<3> >();
    // set index and get centroid coordinate
    grid.setIndex(this->currIdx);
    double xc[3];
    double coordScales[NDIM];
    for (unsigned i = 0; i < NDIM; i++)
    {
      coordScales[i] = 0.5*grid.getDx(i);
    }

    grid.getCentroid(xc);

    // Loop over all node locations on reference element and convert them
    // to appropriate coordinates
    for (unsigned i = 0; i < this->getNumNodes(); i++)
    {
      for (unsigned dim = 0; dim < NDIM; dim++)
      {
        // Node list already has the coordinates of nodes on reference element centered at origin
        nodeCoords(i,dim) = xc[dim] + nodeList(i,dim)*coordScales[dim];
      }
    }
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
    for (unsigned i = 0; i < refNjNk.rows(); i++)
      for (unsigned j = 0; j < refNjNk.cols(); j++)
        NjNk(i,j) = refNjNk(i,j);
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getLowerFaceMassMatrix(unsigned dir,
    Lucee::Matrix<double>& NjNk) const
  {
    if (dir == 0)
    {
      for (unsigned i = 0; i < refFaceNjNk_xl.rows(); i++)
        for (unsigned j = 0; j < refFaceNjNk_xl.cols(); j++)
          NjNk(i,j) = refFaceNjNk_xl(i,j);
    }
    else if (dir == 1)
    {
      for (unsigned i = 0; i < refFaceNjNk_yl.rows(); i++)
        for (unsigned j = 0; j < refFaceNjNk_yl.cols(); j++)
          NjNk(i,j) = refFaceNjNk_yl(i,j);
    }
    else if (dir == 2)
    {
      for (unsigned i = 0; i < refFaceNjNk_zl.rows(); i++)
        for (unsigned j = 0; j < refFaceNjNk_zl.cols(); j++)
          NjNk(i,j) = refFaceNjNk_zl(i,j);
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getUpperFaceMassMatrix(unsigned dir,
    Lucee::Matrix<double>& NjNk) const
  {
    if (dir == 0)
    {
      for (unsigned i = 0; i < refFaceNjNk_xu.rows(); i++)
        for (unsigned j = 0; j < refFaceNjNk_xu.cols(); j++)
          NjNk(i,j) = refFaceNjNk_xu(i,j);
    }
    else if (dir == 1)
    {
      for (unsigned i = 0; i < refFaceNjNk_yu.rows(); i++)
        for (unsigned j = 0; j < refFaceNjNk_yu.cols(); j++)
          NjNk(i,j) = refFaceNjNk_yu(i,j);
    }
    else if (dir == 2)
    {
      for (unsigned i = 0; i < refFaceNjNk_zu.rows(); i++)
        for (unsigned j = 0; j < refFaceNjNk_zu.cols(); j++)
          NjNk(i,j) = refFaceNjNk_zu(i,j);
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getStiffnessMatrix(Lucee::Matrix<double>& DNjDNk) const
  {
    for (unsigned i = 0; i < refDNjDNk.rows(); i++)
        for (unsigned j = 0; j < refDNjDNk.cols(); j++)
          DNjDNk(i,j) = refDNjDNk(i,j);
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::getGradStiffnessMatrix(
    unsigned dir, Lucee::Matrix<double>& DNjNk) const
  {
    if (dir == 0)
    {
      for (unsigned i = 0; i < refDNjNk_0.rows(); i++)
        for (unsigned j = 0; j < refDNjNk_0.cols(); j++)
          DNjNk(i,j) = refDNjNk_0(i,j);
    }
    else if (dir == 1)
    {
      for (unsigned i = 0; i < refDNjNk_1.rows(); i++)
        for (unsigned j = 0; j < refDNjNk_1.cols(); j++)
          DNjNk(i,j) = refDNjNk_1(i,j);
    }
    else if (dir == 2)
    {
      for (unsigned i = 0; i < refDNjNk_2.rows(); i++)
        for (unsigned j = 0; j < refDNjNk_2.cols(); j++)
          DNjNk(i,j) = refDNjNk_2(i,j);
    }
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
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->template getGrid<Lucee::StructuredGridBase<NDIM> >();
    // Get grid spacing (this is assumed to be uniform for now)
    double dx = grid.getDx(0);
    double dy = grid.getDx(1);
    double dz = grid.getDx(2);
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
    // Unused so far: fill out gaussian integration nodes for integration
    // over cell volume
    unsigned totalVolumeGaussNodes = 1.0;
    VectorXd gaussNodeVec(3);
    unsigned currRow = 0;
    double totalWeight;
    /*
    if (NDIM == 3)
    {
      for (unsigned i = 0; i < NDIM; i++)
      {
        totalVolumeGaussNodes = gaussPoints.size()*totalVolumeGaussNodes;
      }

      gaussNodeList = Eigen::MatrixXd(totalVolumeGaussNodes,NDIM+1);
      functionEvaluations  = Eigen::MatrixXd(this->getNumNodes(),totalVolumeGaussNodes);
      
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
            gaussNodeList.row(currRow) << gaussNodeVec(0),gaussNodeVec(1),gaussNodeVec(2),totalWeight;
            for (unsigned functionIndex = 0; functionIndex < functionVector.size(); functionIndex++)
            {
              functionEvaluations(functionIndex,currRow) = evalPolynomial(functionVector[functionIndex],gaussNodeVec);
            }
            currRow++;
          }
        }
      }
    }*/

    
    // Resize the output matrices we need
    resizeMatrices();
    // Call various functions to populate the matrices
    std::cout << "Starting matrix computation" << std::endl;
    computeMass(functionVector,refNjNk);
    computeFaceMass(functionVector,0,refFaceNjNk_xl,refFaceNjNk_xu);
    computeFaceMass(functionVector,1,refFaceNjNk_yl,refFaceNjNk_yu);
    computeFaceMass(functionVector,2,refFaceNjNk_zl,refFaceNjNk_zu);
    computeStiffness(functionVector,refDNjDNk);
    computeGradStiffness(functionVector,0,refDNjNk_0);
    computeGradStiffness(functionVector,1,refDNjNk_1);
    computeGradStiffness(functionVector,2,refDNjNk_2);
    /*std::cout << "Finished computing all matrices" << std::endl;
    std::cout << "refNjNk " << std::endl << refNjNk << std::endl;
    std::cout << "refFaceNjNk_xl " << std::endl << refFaceNjNk_xl << std::endl;
    std::cout << "refFaceNjNk_xu " << std::endl << refFaceNjNk_xu << std::endl;
    std::cout << "refFaceNjNk_yl " << std::endl << refFaceNjNk_yl << std::endl;
    std::cout << "refFaceNjNk_yu " << std::endl << refFaceNjNk_yu << std::endl;
    std::cout << "refFaceNjNk_zl " << std::endl << refFaceNjNk_zl << std::endl;
    std::cout << "refFaceNjNk_zu " << std::endl << refFaceNjNk_zu << std::endl;
    std::cout << "refDNjDNk " << std::endl << refDNjDNk << std::endl;
    std::cout << "refDNjNk_0 " << std::endl << refDNjNk_0 << std::endl;
    std::cout << "refDNjNk_1 "  << std::endl << refDNjNk_1 << std::endl;
    std::cout << "refDNjNk_2 " << std::endl << refDNjNk_2 << std::endl;*/
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
    refDNjNk_0     *= 0.5*dy*0.5*dz;
    refDNjNk_1     *= 0.5*dx*0.5*dz;
    refDNjNk_2     *= 0.5*dx*0.5*dy;
    
    std::cout << "refDNjDNk " << std::endl << refDNjDNk << std::endl;
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::resizeMatrices()
  {
    unsigned generalDim = this->getNumNodes();
    refNjNk        = Eigen::MatrixXd(generalDim,generalDim);
    refDNjDNk      = Eigen::MatrixXd(generalDim,generalDim);
    refDNjNk_0     = Eigen::MatrixXd(generalDim,generalDim);
    refDNjNk_1     = Eigen::MatrixXd(generalDim,generalDim);
    refDNjNk_2     = Eigen::MatrixXd(generalDim,generalDim);
    refFaceNjNk_xl = Eigen::MatrixXd(generalDim,getNumSurfLowerNodes(0));
    refFaceNjNk_xu = Eigen::MatrixXd(generalDim,getNumSurfLowerNodes(0));
    refFaceNjNk_yl = Eigen::MatrixXd(generalDim,getNumSurfLowerNodes(1));
    refFaceNjNk_yu = Eigen::MatrixXd(generalDim,getNumSurfLowerNodes(1));
    refFaceNjNk_zl = Eigen::MatrixXd(generalDim,getNumSurfLowerNodes(2));
    refFaceNjNk_zu = Eigen::MatrixXd(generalDim,getNumSurfLowerNodes(2));
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
  SerendipityElement<NDIM>::evalPolynomial(const blitz::Array<double,3>& polyCoeffs, const VectorXd& nodeCoords)
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
  SerendipityElement<NDIM>::computePolynomialProduct(const blitz::Array<double,3>& poly1, 
    const blitz::Array<double,3>& poly2)
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
  SerendipityElement<NDIM>::computePolynomialDerivative(const blitz::Array<double,3>& poly, unsigned dir)
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
  void
  SerendipityElement<NDIM>::computeMass(const std::vector<blitz::Array<double,3> >& functionVector,
    Eigen::MatrixXd& resultMatrix)
  {
    blitz::Array<double,3> polyProduct(maxPower,maxPower,maxPower);
    blitz::Array<double,2> resultArray(this->getNumNodes(),this->getNumNodes());
    VectorXd gaussNodeVec(3);

    double integrationResult;
    double integrationResultTest;
    double totalWeight;

    resultMatrix.Zero(resultMatrix.rows(),resultMatrix.cols());
    
    for (unsigned kIndex = 0; kIndex < resultMatrix.rows(); kIndex++)
    {
      for (unsigned mIndex = 0; mIndex < resultMatrix.cols(); mIndex++)
      {
        // Reset polyProduct
        polyProduct = 0;
        // Calculate polynomial product of two basis functions
        polyProduct = computePolynomialProduct(functionVector[kIndex],functionVector[mIndex]);

        // Reset integration result
        integrationResult = 0.0;
        integrationResultTest = 0.0;
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
/*
        for (unsigned gaussIndex = 0; gaussIndex < gaussNodeList.rows(); gaussIndex++)
        {
          // result = weight(node) * f(node)getNumSurfLowerNodes()
          integrationResultTest += gaussNodeList(gaussIndex,NDIM+1)*functionEvaluations(kIndex,gaussIndex)*functionEvaluations(mIndex,gaussIndex);
        }

        std::cout << integrationResult - integrationResultTest << std::endl;*/

        resultMatrix(kIndex,mIndex) = integrationResult;
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::computeFaceMass(const std::vector<blitz::Array<double,3> >& functionVector, unsigned dir,
    Eigen::MatrixXd& lowerResultMatrix, Eigen::MatrixXd& upperResultMatrix)
  {
    blitz::Array<double,3> upperPolyProduct(maxPower,maxPower,maxPower);
    blitz::Array<double,3> lowerPolyProduct(maxPower,maxPower,maxPower);
    std::vector<int> surfLowerNodeNums(lowerResultMatrix.cols(),0);
    std::vector<int> surfUpperNodeNums(lowerResultMatrix.cols(),0);
    VectorXd gaussNodeVec(3);

    getSurfLowerNodeNums(dir,surfLowerNodeNums);
    getSurfUpperNodeNums(dir,surfUpperNodeNums);
    
    upperResultMatrix.Zero(upperResultMatrix.rows(),upperResultMatrix.cols());
    lowerResultMatrix.Zero(lowerResultMatrix.rows(),lowerResultMatrix.cols());

    double integrationResultU;
    double integrationResultL;
    double totalWeight;
    unsigned upperNodeNum;
    unsigned lowerNodeNum;

    for (unsigned kIndex = 0; kIndex < upperResultMatrix.rows(); kIndex++)
    {
      // Note that mIndex will be those of the ones in lower/upper node nums
      for (unsigned mIndex = 0; mIndex < upperResultMatrix.cols(); mIndex++)
      {
        // Reset polyProduct
        lowerPolyProduct = 0;
        upperPolyProduct = 0;
        // Assign node nums
        lowerNodeNum = surfLowerNodeNums[mIndex];
        upperNodeNum = surfUpperNodeNums[mIndex];
        // Calculate polynomial product of two basis functions
        lowerPolyProduct = computePolynomialProduct(functionVector[kIndex],functionVector[lowerNodeNum]);
        upperPolyProduct = computePolynomialProduct(functionVector[kIndex],functionVector[upperNodeNum]);

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
            integrationResultU += totalWeight*evalPolynomial(upperPolyProduct,gaussNodeVec);
            // Replace the appropriate coordinate to evaluation location
            gaussNodeVec(dir) = -1;
            integrationResultL += totalWeight*evalPolynomial(lowerPolyProduct,gaussNodeVec);
          }
        }

        upperResultMatrix(kIndex,mIndex) = integrationResultU;
        lowerResultMatrix(kIndex,mIndex) = integrationResultL;
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::computeStiffness(const std::vector<blitz::Array<double,3> >& functionVector,
    Eigen::MatrixXd& resultMatrix)
  {
    blitz::Array<double,3> polyProduct(maxPower,maxPower,maxPower);
    blitz::Array<double,3> functionDx2(maxPower,maxPower,maxPower);
    blitz::Array<double,3> functionDy2(maxPower,maxPower,maxPower);
    blitz::Array<double,3> functionDz2(maxPower,maxPower,maxPower);
    VectorXd gaussNodeVec(3);

    // Get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->template getGrid<Lucee::StructuredGridBase<NDIM> >();
    // Get grid spacing (this is assumed to be uniform for now)
    double dx = grid.getDx(0);
    double dy = grid.getDx(1);
    double dz = grid.getDx(2);
    double dx2 = dx*dx;
    double dy2 = dy*dy;
    double dz2 = dz*dz;

    double integrationResult;
    double totalWeight;

    resultMatrix.Zero(resultMatrix.rows(),resultMatrix.cols());
    
    for (unsigned kIndex = 0; kIndex < resultMatrix.rows(); kIndex++)
    {
      for (unsigned mIndex = 0; mIndex < resultMatrix.cols(); mIndex++)
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
        // Calculate stiffness matrix integrand
        polyProduct = functionDx2*(4.0/dx2) + functionDy2*(4.0/dy2) + functionDz2*(4.0/dz2);

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

        resultMatrix(kIndex,mIndex) = integrationResult;
      }
    }
  }

  template <unsigned NDIM>
  void
  SerendipityElement<NDIM>::computeGradStiffness(const std::vector<blitz::Array<double,3> >& functionVector, 
    unsigned dir, Eigen::MatrixXd& resultMatrix)
  {
    blitz::Array<double,3> polyProduct(maxPower,maxPower,maxPower);
    blitz::Array<double,3> functionD(maxPower,maxPower,maxPower);
    VectorXd gaussNodeVec(3);

    double integrationResult;
    double totalWeight;

    resultMatrix.Zero(resultMatrix.rows(),resultMatrix.cols());

    for (unsigned kIndex = 0; kIndex < resultMatrix.rows(); kIndex++)
    {
      for (unsigned mIndex = 0; mIndex < resultMatrix.cols(); mIndex++)
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

        resultMatrix(kIndex,mIndex) = integrationResult;
      }
    }
  }
  
// instantiations
  template class SerendipityElement<1>;
  template class SerendipityElement<2>;
  template class SerendipityElement<3>;
  template class SerendipityElement<4>;
  template class SerendipityElement<5>;
}