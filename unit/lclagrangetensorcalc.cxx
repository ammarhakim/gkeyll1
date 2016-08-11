/**
 * @file	lclagrangetensorbasiscalc.cxx
 *
 * @brief	Unit tests for Lagrange tensor tests.
 */

// lucee includes
#include <LcLagrangeTensorBasisCalc.h>
#include <LcTest.h>

void
test_0()
{ // default ctor
  Lucee::LagrangeTensorBasisCalc<1> basis;
  LC_ASSERT("Checking number of nodes", basis.getNumNodes() == 1);

  std::vector<double> nloc = basis.getNodeLoc(0);
  LC_ASSERT("Testing node locations", nloc[0] == 0.0);
  
  Lucee::Matrix<double> coeff(1,1);
  basis.getCoeffMat(coeff);
  LC_ASSERT("Testing coefficient matrix", coeff(0,0) == 1.0);
}

void
test_1()
{ // 1D, 2 nodes
  Lucee::LagrangeTensorBasisCalc<1> basis;
  unsigned nn = 2;
  unsigned numNodes[1];
  numNodes[0] = nn;
  basis.calc(Lucee::UNIFORM, numNodes);

  LC_ASSERT("Checking number of nodes", basis.getNumNodes() == nn);

  std::vector<double> nloc = basis.getNodeLoc(0);
  LC_ASSERT("Testing node locations", nloc[0] == -1.0);
  LC_ASSERT("Testing node locations", nloc[1] == 1.0);
  
  Lucee::Matrix<double> coeff(nn,nn);
  basis.getCoeffMat(coeff);
  LC_ASSERT("Testing coefficient matrix", coeff(0,0) == 0.5);
  LC_ASSERT("Testing coefficient matrix", coeff(1,0) == -0.5);

  LC_ASSERT("Testing coefficient matrix", coeff(0,1) == 0.5);
  LC_ASSERT("Testing coefficient matrix", coeff(1,1) == 0.5);

  double xc[1];
  for (unsigned b=0; b<nn; ++b)
  {
    for (unsigned n=0; n<nn; ++n)
    {
      xc[0] = nloc[n];

      if (n == b)
        LC_ASSERT("Testing basis function eval", epsCmp(basis.evalBasis(b, xc), 1.0));
      else
        LC_ASSERT("Testing basis function eval", epsCmp(basis.evalBasis(b, xc), 0.0));
    }
  }

  unsigned tn = basis.getNumNodes();
  Lucee::Matrix<double> massMatrix(tn, tn);
  basis.getMassMatrix(massMatrix);

  LC_ASSERT("Tesing mass-matrix", epsCmp(massMatrix(0,0), 2.0/3.0));
  LC_ASSERT("Tesing mass-matrix", epsCmp(massMatrix(0,1), 1.0/3.0));
  LC_ASSERT("Tesing mass-matrix", epsCmp(massMatrix(1,0), 1.0/3.0));
  LC_ASSERT("Tesing mass-matrix", epsCmp(massMatrix(1,1), 2.0/3.0));

  Lucee::Matrix<double> gradStiffMatrix(tn, tn);
  basis.getGradStiffnessMatrix(0, gradStiffMatrix);

  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(0,0), -1/2.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(0,1), -1/2.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(1,0), 1/2.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(1,1), 1/2.0));

// node lists
  std::vector<int> en;
  basis.getExclusiveNodeIndices(en);

  LC_ASSERT("Testing number of exclusive nodes", en.size() == 1);
  LC_ASSERT("Testing exclusive indices", en[0] == 0);

// list of indices
  std::vector<std::vector<int> > eni;
  basis.getExclusiveNodeNdimIndices(eni);
  
  LC_ASSERT("Testing number of exclusive nodes", eni.size() == 1);
  LC_ASSERT("Testing number of exclusive nodes", eni[0].size() == 1);
  LC_ASSERT("Testing exclusive indices", eni[0][0] == 0);

  std::vector<int> fn;

  LC_ASSERT("Testing number of lower surface nodes", basis.getNumSurfLowerNodes(0) == 1);
  basis.getSurfLowerNodeNums(0, fn);
  LC_ASSERT("Tesing lower node numbers", fn[0] == 0);

  LC_ASSERT("Testing number of upper surface nodes", basis.getNumSurfUpperNodes(0) == 1);
  basis.getSurfUpperNodeNums(0, fn);
  LC_ASSERT("Tesing lower node numbers", fn[0] == 1);

// stiffness matrices
  double dx[4] = {2.0, 2.0, 2.0, 2.0};
  Lucee::Matrix<double> stiff(basis.getNumNodes(), basis.getNumNodes());
  basis.getStiffnessMatrix(dx, stiff);

  LC_ASSERT("Testing stiffness matrix", epsCmp(stiff(0,0), 0.5));
  LC_ASSERT("Testing stiffness matrix", epsCmp(stiff(0,1), -0.5));
  LC_ASSERT("Testing stiffness matrix", epsCmp(stiff(1,0), -0.5));
  LC_ASSERT("Testing stiffness matrix", epsCmp(stiff(1,1), 0.5));

// check quadrature data
  Lucee::Matrix<double> vInterpMat(nn,nn), vOrds(nn,1);
  std::vector<double> vWeights(nn);
  basis.getGaussQuadData(vInterpMat, vOrds, vWeights);

  LC_ASSERT("Testing ordinates", epsCmp(vOrds(0,0), -1/std::sqrt(3)));
  LC_ASSERT("Testing ordinates", epsCmp(vOrds(1,0), 1/std::sqrt(3)));

  LC_ASSERT("Testing weights", epsCmp(vWeights[0], 1.0));
  LC_ASSERT("Testing weights", epsCmp(vWeights[1], 1.0));

  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(0,0), 1/sqrt(3)/2.0+1.0/2.0));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(0,1), 1.0/2.0-1/sqrt(3)/2.0));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(1,0), 1.0/2.0-1/sqrt(3)/2.0));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(1,1), 1/sqrt(3)/2.0+1.0/2.0));

// testing nodal weights
  std::vector<double> nodalWeights(basis.getNumNodes());
  basis.getWeights(nodalWeights);

  LC_ASSERT("Testing nodal weights", epsCmp(nodalWeights[0], 1.0));
  LC_ASSERT("Testing nodal weights", epsCmp(nodalWeights[1], 1.0));
}

void
test_2()
{  // 1D, 3 nodes
  Lucee::LagrangeTensorBasisCalc<1> basis;
  unsigned nn = 3;
  unsigned numNodes[1];
  numNodes[0] = nn;
  basis.calc(Lucee::UNIFORM, numNodes);

  LC_ASSERT("Checking number of nodes", basis.getNumNodes() == nn);

  std::vector<double> nloc = basis.getNodeLoc(0);
  LC_ASSERT("Testing node locations", nloc[0] == -1.0);
  LC_ASSERT("Testing node locations", nloc[1] == 0.0);
  LC_ASSERT("Testing node locations", nloc[2] == 1.0);
  
  Lucee::Matrix<double> coeff(nn,nn);
  basis.getCoeffMat(coeff);

  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(0,0), 1.0/6.0));
  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(1,0), -1.0/2.0));
  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(2,0), 1.0/3.0));

  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(0,1), 2.0/3.0));
  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(1,1), 0.0));
  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(2,1), -2.0/3.0));

  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(0,2), 1.0/6.0));
  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(1,2), 1.0/2.0));
  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(2,2), 1.0/3.0));

  double xc[1];
  for (unsigned b=0; b<nn; ++b)
  {
    for (unsigned n=0; n<nn; ++n)
    {
      xc[0] = nloc[n];

      if (n == b)
        LC_ASSERT("Testing basis function eval", epsCmp(basis.evalBasis(b, xc), 1.0));
      else
        LC_ASSERT("Testing basis function eval", epsCmp(basis.evalBasis(b, xc), 0.0));
    }
  }

  unsigned tn = basis.getNumNodes();
  Lucee::Matrix<double> massMatrix(tn, tn);
  basis.getMassMatrix(massMatrix);

  LC_ASSERT("Tesing mass-matrix", epsCmp(massMatrix(0,0), 4.0/15.0));
  LC_ASSERT("Tesing mass-matrix", epsCmp(massMatrix(0,1), 2.0/15.0));
  LC_ASSERT("Tesing mass-matrix", epsCmp(massMatrix(0,2), (-1.0)/15.0));
  LC_ASSERT("Tesing mass-matrix", epsCmp(massMatrix(1,0), 2.0/15.0));
  LC_ASSERT("Tesing mass-matrix", epsCmp(massMatrix(1,1), 16.0/15.0));
  LC_ASSERT("Tesing mass-matrix", epsCmp(massMatrix(1,2), 2.0/15.0));
  LC_ASSERT("Tesing mass-matrix", epsCmp(massMatrix(2,0), (-1.0)/15.0));
  LC_ASSERT("Tesing mass-matrix", epsCmp(massMatrix(2,1), 2.0/15.0));
  LC_ASSERT("Tesing mass-matrix", epsCmp(massMatrix(2,2), 4.0/15.0));

  Lucee::Matrix<double> gradStiffMatrix(tn, tn);
  basis.getGradStiffnessMatrix(0, gradStiffMatrix);

  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(0,0), -1/2.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(0,1), (-4.0)/(3.0*2.0)));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(0,2), 1/2.0/3.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(1,0), 4.0/(3.0*2.0)));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(1,1), 0.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(1,2), (-4.0)/(3.0*2.0)));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(2,0), -1/2.0/3.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(2,1), 4.0/(3.0*2.0)));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(2,2), 1/2.0));

// node lists
  std::vector<int> en;
  basis.getExclusiveNodeIndices(en);

  LC_ASSERT("Testing number of exclusive nodes", en.size() == 2);
  LC_ASSERT("Testing exclusive indices", en[0] == 0);
  LC_ASSERT("Testing exclusive indices", en[1] == 1);

// list of indices
  std::vector<std::vector<int> > eni;
  basis.getExclusiveNodeNdimIndices(eni);
  
  LC_ASSERT("Testing number of exclusive nodes", eni.size() == 2);
  LC_ASSERT("Testing number of exclusive nodes", eni[0].size() == 1);
  LC_ASSERT("Testing exclusive indices", eni[0][0] == 0);
  LC_ASSERT("Testing exclusive indices", eni[1][0] == 1);

  std::vector<int> fn;

  LC_ASSERT("Testing number of lower surface nodes", basis.getNumSurfLowerNodes(0) == 1);
  basis.getSurfLowerNodeNums(0, fn);
  LC_ASSERT("Tesing lower node numbers", fn[0] == 0);

  LC_ASSERT("Testing number of upper surface nodes", basis.getNumSurfUpperNodes(0) == 1);
  basis.getSurfUpperNodeNums(0, fn);
  LC_ASSERT("Tesing lower node numbers", fn[0] == 2);

// stiffness matrices
  double dx[4] = {2.0, 2.0, 2.0, 2.0};
  Lucee::Matrix<double> stiff(basis.getNumNodes(), basis.getNumNodes());
  basis.getStiffnessMatrix(dx, stiff);

  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(0,0), 7.0/6.0, 2e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(0,1), (-4.0)/3.0, 2e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(0,2), 1.0/6.0, 2e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(1,0), (-4.0)/3.0, 2e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(1,1), 8.0/3.0, 2e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(1,2), (-4.0)/3.0, 2e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(2,0), 1.0/6.0, 2e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(2,1), (-4.0)/3.0, 2e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(2,2), 7.0/6.0, 2e-15));

// check quadrature data
  Lucee::Matrix<double> vInterpMat(nn,nn), vOrds(nn,1);
  std::vector<double> vWeights(nn);
  basis.getGaussQuadData(vInterpMat, vOrds, vWeights);

  LC_ASSERT("Testing ordinates", epsCmp(vOrds(0,0), -std::sqrt(3)/std::sqrt(5)));
  LC_ASSERT("Testing ordinates", epsCmp(vOrds(1,0), 0.0));
  LC_ASSERT("Testing ordinates", epsCmp(vOrds(2,0),  std::sqrt(3)/std::sqrt(5)));

  LC_ASSERT("Testing weights", epsCmp(vWeights[0], 5.0/9.0));
  LC_ASSERT("Testing weights", epsCmp(vWeights[1], 8.0/9.0));
  LC_ASSERT("Testing weights", epsCmp(vWeights[2], 5.0/9.0));

  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(0,0), .6872983346207417));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(0,1), .3999999999999997));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(0,2), -.08729833462074155, 10));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(1,0), 0.0));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(1,1), 1.0));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(1,2), 0.0));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(2,0), -.08729833462074174, 10));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(2,1), .4000000000000001));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(2,2), .6872983346207415));

// testing nodal weights
  std::vector<double> nodalWeights(basis.getNumNodes());
  basis.getWeights(nodalWeights);

  LC_ASSERT("Testing nodal weights", epsCmp(nodalWeights[0], 1.0/3.0));
  LC_ASSERT("Testing nodal weights", epsCmp(nodalWeights[1], 4.0/3.0));
  LC_ASSERT("Testing nodal weights", epsCmp(nodalWeights[2], 1.0/3.0));
}

void
test_2_lobatto()
{  // 1D, 3 nodes
  Lucee::LagrangeTensorBasisCalc<1> basis;
  unsigned nn = 3;
  unsigned numNodes[1];
  numNodes[0] = nn;
  basis.calc(Lucee::LOBATTO, numNodes);

  LC_ASSERT("Checking number of nodes", basis.getNumNodes() == nn);

  std::vector<double> nloc = basis.getNodeLoc(0);
  LC_ASSERT("Testing node locations", nloc[0] == -1.0);
  LC_ASSERT("Testing node locations", nloc[1] == 0.0);
  LC_ASSERT("Testing node locations", nloc[2] == 1.0);
  
  Lucee::Matrix<double> coeff(nn,nn);
  basis.getCoeffMat(coeff);

  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(0,0), 1.0/6.0));
  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(1,0), -1.0/2.0));
  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(2,0), 1.0/3.0));

  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(0,1), 2.0/3.0));
  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(1,1), 0.0));
  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(2,1), -2.0/3.0));

  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(0,2), 1.0/6.0));
  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(1,2), 1.0/2.0));
  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(2,2), 1.0/3.0));

  double xc[1];
  for (unsigned b=0; b<nn; ++b)
  {
    for (unsigned n=0; n<nn; ++n)
    {
      xc[0] = nloc[n];

      if (n == b)
        LC_ASSERT("Testing basis function eval", epsCmp(basis.evalBasis(b, xc), 1.0));
      else
        LC_ASSERT("Testing basis function eval", epsCmp(basis.evalBasis(b, xc), 0.0));
    }
  }

  unsigned tn = basis.getNumNodes();
  Lucee::Matrix<double> massMatrix(tn, tn);
  basis.getMassMatrix(massMatrix);

  LC_ASSERT("Tesing mass-matrix", epsCmp(massMatrix(0,0), 4.0/15.0));
  LC_ASSERT("Tesing mass-matrix", epsCmp(massMatrix(0,1), 2.0/15.0));
  LC_ASSERT("Tesing mass-matrix", epsCmp(massMatrix(0,2), (-1.0)/15.0));
  LC_ASSERT("Tesing mass-matrix", epsCmp(massMatrix(1,0), 2.0/15.0));
  LC_ASSERT("Tesing mass-matrix", epsCmp(massMatrix(1,1), 16.0/15.0));
  LC_ASSERT("Tesing mass-matrix", epsCmp(massMatrix(1,2), 2.0/15.0));
  LC_ASSERT("Tesing mass-matrix", epsCmp(massMatrix(2,0), (-1.0)/15.0));
  LC_ASSERT("Tesing mass-matrix", epsCmp(massMatrix(2,1), 2.0/15.0));
  LC_ASSERT("Tesing mass-matrix", epsCmp(massMatrix(2,2), 4.0/15.0));

  Lucee::Matrix<double> gradStiffMatrix(tn, tn);
  basis.getGradStiffnessMatrix(0, gradStiffMatrix);

  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(0,0), -1/2.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(0,1), (-4.0)/(3.0*2.0)));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(0,2), 1/2.0/3.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(1,0), 4.0/(3.0*2.0)));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(1,1), 0.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(1,2), (-4.0)/(3.0*2.0)));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(2,0), -1/2.0/3.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(2,1), 4.0/(3.0*2.0)));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(2,2), 1/2.0));

// node lists
  std::vector<int> en;
  basis.getExclusiveNodeIndices(en);

  LC_ASSERT("Testing number of exclusive nodes", en.size() == 2);
  LC_ASSERT("Testing exclusive indices", en[0] == 0);
  LC_ASSERT("Testing exclusive indices", en[1] == 1);

  std::vector<int> fn;

  LC_ASSERT("Testing number of lower surface nodes", basis.getNumSurfLowerNodes(0) == 1);
  basis.getSurfLowerNodeNums(0, fn);
  LC_ASSERT("Tesing lower node numbers", fn[0] == 0);

  LC_ASSERT("Testing number of upper surface nodes", basis.getNumSurfUpperNodes(0) == 1);
  basis.getSurfUpperNodeNums(0, fn);
  LC_ASSERT("Tesing lower node numbers", fn[0] == 2);

// check quadrature data
  Lucee::Matrix<double> vInterpMat(nn,nn), vOrds(nn,1);
  std::vector<double> vWeights(nn);
  basis.getGaussQuadData(vInterpMat, vOrds, vWeights);

  LC_ASSERT("Testing ordinates", epsCmp(vOrds(0,0), -std::sqrt(3)/std::sqrt(5)));
  LC_ASSERT("Testing ordinates", epsCmp(vOrds(1,0), 0.0));
  LC_ASSERT("Testing ordinates", epsCmp(vOrds(2,0),  std::sqrt(3)/std::sqrt(5)));

  LC_ASSERT("Testing weights", epsCmp(vWeights[0], 5.0/9.0));
  LC_ASSERT("Testing weights", epsCmp(vWeights[1], 8.0/9.0));
  LC_ASSERT("Testing weights", epsCmp(vWeights[2], 5.0/9.0));

  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(0,0), .6872983346207417));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(0,1), .3999999999999997));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(0,2), -.08729833462074155, 10));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(1,0), 0.0));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(1,1), 1.0));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(1,2), 0.0));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(2,0), -.08729833462074174, 10));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(2,1), .4000000000000001));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(2,2), .6872983346207415));
}

void
test_3()
{  // 1D, 4 nodes
  Lucee::LagrangeTensorBasisCalc<1> basis;
  unsigned nn = 4;
  unsigned numNodes[1];
  numNodes[0] = nn;
  basis.calc(Lucee::UNIFORM, numNodes);

  LC_ASSERT("Checking number of nodes", basis.getNumNodes() == nn);

  std::vector<double> nloc = basis.getNodeLoc(0);
  LC_ASSERT("Testing node locations", nloc[0] == -1.0);
  LC_ASSERT("Testing node locations", epsCmp(nloc[1], -1.0/3.0));
  LC_ASSERT("Testing node locations", epsCmp(nloc[2], 1.0/3.0));
  LC_ASSERT("Testing node locations", nloc[3] == 1.0);
  
  Lucee::Matrix<double> coeff(nn,nn);
  basis.getCoeffMat(coeff);

  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(0,0), 1./8));
  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(1,0), -11./40));
  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(2,0), 3./8));
  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(3,0), -9./40));

  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(0,1), 3./8));
  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(1,1), -27./40));
  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(2,1), -3./8));
  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(3,1), 27./40));

  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(0,2), 3./8));
  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(1,2), 27./40));
  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(2,2), -3./8));
  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(3,2), -27./40));

  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(0,3), 1./8));
  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(1,3), 11./40));
  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(2,3), 3./8));
  LC_ASSERT("Testing coefficient matrix", epsCmp(coeff(3,3), 9./40));

  double xc[1];
  for (unsigned b=0; b<nn; ++b)
  {
    for (unsigned n=0; n<nn; ++n)
    {
      xc[0] = nloc[n];

      if (n == b)
        LC_ASSERT("Testing basis function eval", epsCmp(basis.evalBasis(b, xc), 1.0));
      else
        LC_ASSERT("Testing basis function eval", epsCmp(basis.evalBasis(b, xc), 0.0));
    }
  }

  unsigned tn = basis.getNumNodes();
  Lucee::Matrix<double> massMatrix(tn, tn);
  basis.getMassMatrix(massMatrix);

// NOTE of 10/17/2012: It is possible that some tests below fail due
// to precision level errors. For some reason a few values are
// slightly off (but with relative error of 1e-13 or so), however this
// should not be a problem. (Ammar Hakim)

  LC_ASSERT("Testing mass matrix", epsCmp(massMatrix(0,0), 16.0/105.0));
  LC_ASSERT("Testing mass matrix", epsCmp(massMatrix(0,1), 33.0/280.0));
  LC_ASSERT("Testing mass matrix", epsCmp(massMatrix(0,2), (-3.0)/70.0, 10));
  LC_ASSERT("Testing mass matrix", epsCmp(massMatrix(0,3), 19.0/840.0));
  LC_ASSERT("Testing mass matrix", epsCmp(massMatrix(1,0), 33.0/280.0));
  LC_ASSERT("Testing mass matrix", epsCmp(massMatrix(1,1), 27.0/35.0));
  LC_ASSERT("Testing mass matrix", epsCmp(massMatrix(1,2), (-27.0)/280.0));
  LC_ASSERT("Testing mass matrix", epsCmp(massMatrix(1,3), (-3.0)/70.0));
  LC_ASSERT("Testing mass matrix", epsCmp(massMatrix(2,0), (-3.0)/70.0, 10));
  LC_ASSERT("Testing mass matrix", epsCmp(massMatrix(2,1), (-27.0)/280.0));
  LC_ASSERT("Testing mass matrix", epsCmp(massMatrix(2,2), 27.0/35.0));
  LC_ASSERT("Testing mass matrix", epsCmp(massMatrix(2,3), 33.0/280.0));
  LC_ASSERT("Testing mass matrix", epsCmp(massMatrix(3,0), 19.0/840.0));
  LC_ASSERT("Testing mass matrix", epsCmp(massMatrix(3,1), (-3.0)/70.0));
  LC_ASSERT("Testing mass matrix", epsCmp(massMatrix(3,2), 33.0/280.0));
  LC_ASSERT("Testing mass matrix", epsCmp(massMatrix(3,3), 16.0/105.0));

  Lucee::Matrix<double> gradStiffMatrix(tn, tn);
  basis.getGradStiffnessMatrix(0, gradStiffMatrix);

  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(0,0), -1/2.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(0,1), (-57.0)/(40.0*2.0)));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(0,2), 3.0/(5.0*2.0)));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(0,3), (-7.0)/(40.0*2.0), 10));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(1,0), 57.0/(40.0*2.0)));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(1,1), 0.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(1,2), (-81.0)/(40.0*2.0)));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(1,3), 3.0/(5.0*2.0)));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(2,0), (-3.0)/(5.0*2.0)));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(2,1), 81.0/(40.0*2.0)));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(2,2), 0.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(2,3), (-57.0)/(40.0*2.0)));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(3,0), 7.0/(40.0*2.0)));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(3,1), (-3.0)/(5.0*2.0)));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(3,2), 57.0/(40.0*2.0)));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(3,3), 1/2.0));

// node lists
  std::vector<int> en;
  basis.getExclusiveNodeIndices(en);

  LC_ASSERT("Testing number of exclusive nodes", en.size() == 3);
  LC_ASSERT("Testing exclusive indices", en[0] == 0);
  LC_ASSERT("Testing exclusive indices", en[1] == 1);
  LC_ASSERT("Testing exclusive indices", en[2] == 2);

// list of indices
  std::vector<std::vector<int> > eni;
  basis.getExclusiveNodeNdimIndices(eni);
  
  LC_ASSERT("Testing number of exclusive nodes", eni.size() == 3);
  LC_ASSERT("Testing exclusive indices", eni[0][0] == 0);
  LC_ASSERT("Testing exclusive indices", eni[1][0] == 1);
  LC_ASSERT("Testing exclusive indices", eni[2][0] == 2);

  std::vector<int> fn;

  LC_ASSERT("Testing number of lower surface nodes", basis.getNumSurfLowerNodes(0) == 1);
  basis.getSurfLowerNodeNums(0, fn);
  LC_ASSERT("Tesing lower node numbers", fn[0] == 0);

  LC_ASSERT("Testing number of upper surface nodes", basis.getNumSurfUpperNodes(0) == 1);
  basis.getSurfUpperNodeNums(0, fn);
  LC_ASSERT("Tesing lower node numbers", fn[0] == 3);

// stiffness matrices
  double dx[4] = {2.0, 2.0, 2.0, 2.0};
  Lucee::Matrix<double> stiff(basis.getNumNodes(), basis.getNumNodes());
  basis.getStiffnessMatrix(dx, stiff);

  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(0,0), 37.0/20.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(0,1), (-189.0)/80.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(0,2), 27.0/40.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(0,3), (-13.0)/80.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(1,0), (-189.0)/80.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(1,1), 27.0/5.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(1,2), (-297.0)/80.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(1,3), 27.0/40.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(2,0), 27.0/40.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(2,1), (-297.0)/80.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(2,2), 27.0/5.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(2,3), (-189.0)/80.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(3,0), (-13.0)/80.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(3,1), 27.0/40.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(3,2), (-189.0)/80.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(3,3), 37.0/20.0, 5e-15));

// check quadrature data
  Lucee::Matrix<double> vInterpMat(nn,nn), vOrds(nn,1);
  std::vector<double> vWeights(nn);
  basis.getGaussQuadData(vInterpMat, vOrds, vWeights);

  double x1 = std::sqrt((3-2*sqrt(6./5.))/7);
  double x2 = std::sqrt((3+2*sqrt(6./5.))/7);

  LC_ASSERT("Testing ordinates", epsCmp(vOrds(0,0), -x2));
  LC_ASSERT("Testing ordinates", epsCmp(vOrds(1,0), -x1));
  LC_ASSERT("Testing ordinates", epsCmp(vOrds(2,0), x1));
  LC_ASSERT("Testing ordinates", epsCmp(vOrds(3,0), x2));

  LC_ASSERT("Testing weights", epsCmp(vWeights[0], (18-sqrt(30))/36));
  LC_ASSERT("Testing weights", epsCmp(vWeights[1], (18+sqrt(30))/36));
  LC_ASSERT("Testing weights", epsCmp(vWeights[2], (18+sqrt(30))/36));
  LC_ASSERT("Testing weights", epsCmp(vWeights[3], (18-sqrt(30))/36));

  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(0,0), .6600056650728031, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(0,1), .5209376877117045, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(0,2), -.2301879032507395, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(0,3), 0.0492445504662319, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(1,0), .003373736432772362, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(1,1), 1.004885854825647, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(1,2), -.009921353572325708, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(1,3), .001661762313906814, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(2,0), .001661762313906398, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(2,1), -0.00992135357232471, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(2,2), 1.004885854825646, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(2,3), .003373736432772445, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(3,0), .04924455046623186, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(3,1), -.2301879032507392, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(3,2), 0.520937687711704, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(3,3), .6600056650728033, 2e-15));
}

void
test_4()
{  // 1D, 4 nodes
  Lucee::LagrangeTensorBasisCalc<1> basis;
  unsigned nn = 4;
  unsigned numNodes[1];
  numNodes[0] = nn;
  basis.calc(Lucee::GAUSSIAN, numNodes);

  LC_ASSERT("Checking number of nodes", basis.getNumNodes() == nn);

  double x1 = std::sqrt((3-2*sqrt(6./5.))/7);
  double x2 = std::sqrt((3+2*sqrt(6./5.))/7);

  std::vector<double> nloc = basis.getNodeLoc(0);
  LC_ASSERT("Testing node locations", epsCmp(nloc[0], -x2));
  LC_ASSERT("Testing node locations", epsCmp(nloc[1], -x1));
  LC_ASSERT("Testing node locations", epsCmp(nloc[2], x1));
  LC_ASSERT("Testing node locations", epsCmp(nloc[3], x2));

// I do not have an explicit test for the coefficient matrix. (Too
// lazy to compute them by hand). However, it is tested via the
// evalBasis function below. (Ammar Hakim, 10/16/2012)
  
  double xc[1];
  for (unsigned b=0; b<nn; ++b)
  {
    for (unsigned n=0; n<nn; ++n)
    {
      xc[0] = nloc[n];

      if (n == b)
        LC_ASSERT("Testing basis function eval", epsCmp(basis.evalBasis(b, xc), 1.0));
      else
        LC_ASSERT("Testing basis function eval", epsCmp(basis.evalBasis(b, xc), 0.0));
    }
  }

// node lists
  std::vector<int> en;
  basis.getExclusiveNodeIndices(en);

  LC_ASSERT("Testing number of exclusive nodes", en.size() == 4);
  LC_ASSERT("Testing exclusive indices", en[0] == 0);
  LC_ASSERT("Testing exclusive indices", en[1] == 1);
  LC_ASSERT("Testing exclusive indices", en[2] == 2);
  LC_ASSERT("Testing exclusive indices", en[3] == 3);

  LC_ASSERT("Testing number of lower surface nodes", basis.getNumSurfLowerNodes(0) == 0);
  LC_ASSERT("Testing number of upper surface nodes", basis.getNumSurfUpperNodes(0) == 0);

// check quadrature data
  Lucee::Matrix<double> vInterpMat(nn,nn), vOrds(nn,1);
  std::vector<double> vWeights(nn);
  basis.getGaussQuadData(vInterpMat, vOrds, vWeights);

  LC_ASSERT("Testing ordinates", epsCmp(vOrds(0,0), -x2));
  LC_ASSERT("Testing ordinates", epsCmp(vOrds(1,0), -x1));
  LC_ASSERT("Testing ordinates", epsCmp(vOrds(2,0), x1));
  LC_ASSERT("Testing ordinates", epsCmp(vOrds(3,0), x2));

  LC_ASSERT("Testing weights", epsCmp(vWeights[0], (18-sqrt(30))/36));
  LC_ASSERT("Testing weights", epsCmp(vWeights[1], (18+sqrt(30))/36));
  LC_ASSERT("Testing weights", epsCmp(vWeights[2], (18+sqrt(30))/36));
  LC_ASSERT("Testing weights", epsCmp(vWeights[3], (18-sqrt(30))/36));

  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(0,0), 1.0));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(0,1), 0.0));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(0,2), 0.0));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(0,3), 0.0));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(1,0), 0.0));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(1,1), 1.0));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(1,2), 0.0));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(1,3), 0.0));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(2,0), 0.0));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(2,1), 0.0));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(2,2), 1.0));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(2,3), 0.0));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(3,0), 0.0));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(3,1), 0.0));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(3,2), 0.0));
  LC_ASSERT("Testing interpolation matrix", epsCmp(vInterpMat(3,3), 1.0));

// testing nodal weights
  std::vector<double> nodalWeights(basis.getNumNodes());
  basis.getWeights(nodalWeights);

  LC_ASSERT("Testing nodal weights", epsCmp(nodalWeights[0], (18-sqrt(30))/36));
  LC_ASSERT("Testing nodal weights", epsCmp(nodalWeights[1], (18+sqrt(30))/36));
  LC_ASSERT("Testing nodal weights", epsCmp(nodalWeights[2], (18+sqrt(30))/36));
  LC_ASSERT("Testing nodal weights", epsCmp(nodalWeights[3], (18-sqrt(30))/36));
}

void
test_4_lobatto()
{  // 1D, 4 nodes
  Lucee::LagrangeTensorBasisCalc<1> basis;
  unsigned nn = 4;
  unsigned numNodes[1];
  numNodes[0] = nn;
  basis.calc(Lucee::LOBATTO, numNodes);

  LC_ASSERT("Checking number of nodes", basis.getNumNodes() == nn);

  double x1 = std::sqrt(1.0/5.0);
  double x2 = 1.0;

  std::vector<double> nloc = basis.getNodeLoc(0);
  LC_ASSERT("Testing node locations", epsCmp(nloc[0], -x2));
  LC_ASSERT("Testing node locations", epsCmp(nloc[1], -x1));
  LC_ASSERT("Testing node locations", epsCmp(nloc[2], x1));
  LC_ASSERT("Testing node locations", epsCmp(nloc[3], x2));

// I do not have an explicit test for the coefficient matrix. (Too
// lazy to compute them by hand). However, it is tested via the
// evalBasis function below. (Ammar Hakim, 10/16/2012)
  
  double xc[1];
  for (unsigned b=0; b<nn; ++b)
  {
    for (unsigned n=0; n<nn; ++n)
    {
      xc[0] = nloc[n];

      if (n == b)
        LC_ASSERT("Testing basis function eval", epsCmp(basis.evalBasis(b, xc), 1.0));
      else
        LC_ASSERT("Testing basis function eval", epsCmp(basis.evalBasis(b, xc), 0.0));
    }
  }

// node lists
  std::vector<int> en;
  basis.getExclusiveNodeIndices(en);

  LC_ASSERT("Testing number of exclusive nodes", en.size() == 3);
  LC_ASSERT("Testing exclusive indices", en[0] == 0);
  LC_ASSERT("Testing exclusive indices", en[1] == 1);
  LC_ASSERT("Testing exclusive indices", en[2] == 2);

  LC_ASSERT("Testing number of lower surface nodes", basis.getNumSurfLowerNodes(0) == 1);
  LC_ASSERT("Testing number of upper surface nodes", basis.getNumSurfUpperNodes(0) == 1);

// check quadrature data
  Lucee::Matrix<double> vInterpMat(nn,nn), vOrds(nn,1);
  std::vector<double> vWeights(nn);
  basis.getGaussQuadData(vInterpMat, vOrds, vWeights);

  x1 = std::sqrt((3-2*sqrt(6./5.))/7);
  x2 = std::sqrt((3+2*sqrt(6./5.))/7);

  LC_ASSERT("Testing ordinates", epsCmp(vOrds(0,0), -x2));
  LC_ASSERT("Testing ordinates", epsCmp(vOrds(1,0), -x1));
  LC_ASSERT("Testing ordinates", epsCmp(vOrds(2,0), x1));
  LC_ASSERT("Testing ordinates", epsCmp(vOrds(3,0), x2));

  LC_ASSERT("Testing weights", epsCmp(vWeights[0], (18-sqrt(30))/36));
  LC_ASSERT("Testing weights", epsCmp(vWeights[1], (18+sqrt(30))/36));
  LC_ASSERT("Testing weights", epsCmp(vWeights[2], (18+sqrt(30))/36));
  LC_ASSERT("Testing weights", epsCmp(vWeights[3], (18-sqrt(30))/36));

// testing nodal weights
  std::vector<double> nodalWeights(basis.getNumNodes());
  basis.getWeights(nodalWeights);

  LC_ASSERT("Testing nodal weights", epsCmp(nodalWeights[0], .1666666666666667));
  LC_ASSERT("Testing nodal weights", epsCmp(nodalWeights[1], .8333333333333333));
  LC_ASSERT("Testing nodal weights", epsCmp(nodalWeights[2], .8333333333333333));
  LC_ASSERT("Testing nodal weights", epsCmp(nodalWeights[3], .1666666666666667));
}

void
test_5()
{ // 2D, 2x2 nodes
  Lucee::LagrangeTensorBasisCalc<2> basis;
  unsigned nn = 2;
  unsigned numNodes[2];
  numNodes[0] = nn; numNodes[1] = nn;
  basis.calc(Lucee::UNIFORM, numNodes);

  LC_ASSERT("Checking number of nodes", basis.getNumNodes() == nn*nn);
  
  double xn[2];
  double loc1[2] = {-1, 1};

  unsigned nIdx = 0;
  for (unsigned i=0; i<nn; ++i)
  {
    for (unsigned j=0; j<nn; ++j)
    {
      basis.fillWithNodeCoordinate(nIdx, xn);

      LC_ASSERT("Testing node coordinate", epsCmp(xn[0], loc1[i]));
      LC_ASSERT("Testing node coordinate", epsCmp(xn[1], loc1[j]));

      nIdx += 1;
    }
  }

  nIdx = 0;
  for (unsigned i=0; i<nn; ++i)
  {
    for (unsigned j=0; j<nn; ++j)
    {
      basis.fillWithNodeCoordinate(nIdx, xn);

      for (unsigned b=0; b<nn*nn; ++b)
      {
        if (nIdx == b)
          LC_ASSERT("Testing basis function eval", epsCmp(basis.evalBasis(b, xn), 1.0));
        else
          LC_ASSERT("Testing basis function eval", epsCmp(basis.evalBasis(b, xn), 0.0));
      }

      nIdx += 1;
    }
  }

  unsigned tn = basis.getNumNodes();
  Lucee::Matrix<double> massMatrix(tn, tn);
  basis.getMassMatrix(massMatrix);

  LC_ASSERT("Testing mass-matrix", epsCmp(massMatrix(0,0), 4.0/9.0));
  LC_ASSERT("Testing mass-matrix", epsCmp(massMatrix(0,1), 2.0/9.0));
  LC_ASSERT("Testing mass-matrix", epsCmp(massMatrix(0,2), 2.0/9.0));
  LC_ASSERT("Testing mass-matrix", epsCmp(massMatrix(0,3), 1.0/9.0));

  LC_ASSERT("Testing mass-matrix", epsCmp(massMatrix(1,0), 2.0/9.0));
  LC_ASSERT("Testing mass-matrix", epsCmp(massMatrix(1,1), 4.0/9.0));
  LC_ASSERT("Testing mass-matrix", epsCmp(massMatrix(1,2), 1.0/9.0));
  LC_ASSERT("Testing mass-matrix", epsCmp(massMatrix(1,3), 2.0/9.0));

  LC_ASSERT("Testing mass-matrix", epsCmp(massMatrix(2,0), 2.0/9.0));
  LC_ASSERT("Testing mass-matrix", epsCmp(massMatrix(2,1), 1.0/9.0));
  LC_ASSERT("Testing mass-matrix", epsCmp(massMatrix(2,2), 4.0/9.0));
  LC_ASSERT("Testing mass-matrix", epsCmp(massMatrix(2,3), 2.0/9.0));

  LC_ASSERT("Testing mass-matrix", epsCmp(massMatrix(3,0), 1.0/9.0));
  LC_ASSERT("Testing mass-matrix", epsCmp(massMatrix(3,1), 2.0/9.0));
  LC_ASSERT("Testing mass-matrix", epsCmp(massMatrix(3,2), 2.0/9.0));
  LC_ASSERT("Testing mass-matrix", epsCmp(massMatrix(3,3), 4.0/9.0));

  Lucee::Matrix<double> gradStiffMatrix(tn, tn);

  basis.getGradStiffnessMatrix(0, gradStiffMatrix);
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(0,0), (-1.0)/3.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(0,1), (-1.0)/6.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(0,2), (-1.0)/3.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(0,3), (-1.0)/6.0));

  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(1,0), (-1.0)/6.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(1,1), (-1.0)/3.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(1,2), (-1.0)/6.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(1,3), (-1.0)/3.0));

  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(2,0), 1.0/3.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(2,1), 1.0/6.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(2,2), 1.0/3.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(2,3), 1.0/6.0));

  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(3,0), 1.0/6.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(3,1), 1.0/3.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(3,2), 1.0/6.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(3,3), 1.0/3.0));

  basis.getGradStiffnessMatrix(1, gradStiffMatrix);
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(0,0), (-2.0)/(3.0*2.0)));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(0,1), (-2.0)/(3.0*2.0)));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(0,2), -1/2.0/3.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(0,3), -1/2.0/3.0));

  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(1,0), 2.0/(3.0*2.0)));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(1,1), 2.0/(3.0*2.0)));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(1,2), 1/2.0/3.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(1,3), 1/2.0/3.0));

  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(2,0), -1/2.0/3.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(2,1), -1/2.0/3.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(2,2), (-2.0)/(3.0*2.0)));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(2,3), (-2.0)/(3.0*2.0)));

  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(3,0), 1/2.0/3.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(3,1), 1/2.0/3.0));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(3,2), 2.0/(3.0*2.0)));
  LC_ASSERT("Testing grad-stiffness matrix", epsCmp(gradStiffMatrix(3,3), 2.0/(3.0*2.0)));

// node lists
  std::vector<int> en;
  basis.getExclusiveNodeIndices(en);

  LC_ASSERT("Testing number of exclusive nodes", en.size() == 1);
  LC_ASSERT("Testing exclusive indices", en[0] == 0);

// list of indices
  std::vector<std::vector<int> > eni;
  basis.getExclusiveNodeNdimIndices(eni);
  
  LC_ASSERT("Testing number of exclusive nodes", eni.size() == 1);
  LC_ASSERT("Testing number of exclusive nodes", eni[0].size() == 2);
  LC_ASSERT("Testing exclusive indices", eni[0][0] == 0);
  LC_ASSERT("Testing exclusive indices", eni[0][1] == 0);

  std::vector<int> fn;

  LC_ASSERT("Testing number of lower surface nodes", basis.getNumSurfLowerNodes(0) == 2);
  basis.getSurfLowerNodeNums(0, fn);
  LC_ASSERT("Tesing lower node numbers", fn[0] == 0);
  LC_ASSERT("Tesing lower node numbers", fn[1] == 1);

  LC_ASSERT("Testing number of lower surface nodes", basis.getNumSurfLowerNodes(1) == 2);
  basis.getSurfLowerNodeNums(1, fn);
  LC_ASSERT("Tesing lower node numbers", fn[0] == 0);
  LC_ASSERT("Tesing lower node numbers", fn[1] == 2);

  LC_ASSERT("Testing number of upper surface nodes", basis.getNumSurfUpperNodes(0) == 2);
  basis.getSurfUpperNodeNums(0, fn);
  LC_ASSERT("Tesing lower node numbers", fn[0] == 2);
  LC_ASSERT("Tesing lower node numbers", fn[1] == 3);

  LC_ASSERT("Testing number of upper surface nodes", basis.getNumSurfUpperNodes(1) == 2);
  basis.getSurfUpperNodeNums(1, fn);
  LC_ASSERT("Tesing lower node numbers", fn[0] == 1);
  LC_ASSERT("Tesing lower node numbers", fn[1] == 3);

// test face mass matrices
  Lucee::Matrix<double> faceMass(basis.getNumNodes(), basis.getNumSurfLowerNodes(0));

  basis.getLowerFaceMassMatrix(0, faceMass);
  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(0,0), 2.0/3.0));
  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(0,1), 1.0/3.0));

  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(1,0), 1.0/3.0));
  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(1,1), 2.0/3.0));

  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(2,0), 0.0));
  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(2,1), 0.0));

  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(3,0), 0.0));
  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(3,1), 0.0));

  basis.getUpperFaceMassMatrix(0, faceMass);
  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(0,0), 0.0));
  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(0,1), 0.0));

  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(1,0), 0.0));
  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(1,1), 0.0));

  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(2,0), 2.0/3.0));
  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(2,1), 1.0/3.0));

  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(3,0), 1.0/3.0));
  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(3,1), 2.0/3.0));

  basis.getLowerFaceMassMatrix(1, faceMass);
  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(0,0), 2.0/3.0));
  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(0,1), 1.0/3.0));

  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(1,0), 0.0));
  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(1,1), 0.0));

  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(2,0), 1.0/3.0));
  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(2,1), 2.0/3.0));

  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(3,0), 0.0));
  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(3,1), 0.0));

  basis.getUpperFaceMassMatrix(1, faceMass);
  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(0,0), 0.0));
  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(0,1), 0.0));

  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(1,0), 2.0/3.0));
  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(1,1), 1.0/3.0));

  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(2,0), 0.0));
  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(2,1), 0.0));

  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(3,0), 1.0/3.0));
  LC_ASSERT("Testing face-mass matrix", epsCmp(faceMass(3,1), 2.0/3.0));

// stiffness matrices
  double dx[4] = {2.0, 2.0, 2.0, 2.0};
  Lucee::Matrix<double> stiff(basis.getNumNodes(), basis.getNumNodes());
  basis.getStiffnessMatrix(dx, stiff);
 
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(0,0), 2.0/3.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(0,1), (-1.0)/6.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(0,2), (-1.0)/6.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(0,3), (-1.0)/3.0, 5e-15));

  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(1,0), (-1.0)/6.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(1,1), 2.0/3.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(1,2), (-1.0)/3.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(1,3), (-1.0)/6.0, 5e-15));

  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(2,0), (-1.0)/6.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(2,1), (-1.0)/3.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(2,2), 2.0/3.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(2,3), (-1.0)/6.0, 5e-15));

  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(3,0), (-1.0)/3.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(3,1), (-1.0)/6.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(3,2), (-1.0)/6.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(3,3), 2.0/3.0, 5e-15));

  unsigned numTotalNodes = basis.getNumNodes();
// check quadrature data
  Lucee::Matrix<double> vInterpMat(numTotalNodes,numTotalNodes), vOrds(numTotalNodes,2);
  std::vector<double> vWeights(numTotalNodes);
  basis.getGaussQuadData(vInterpMat, vOrds, vWeights);

  double ords[2] = {-1/std::sqrt(3), 1/std::sqrt(3)};
  double weights[2] = {1.0, 1.0};

  LC_ASSERT("Testing ordinates", epsCmp(vOrds(0,0), ords[0]));
  LC_ASSERT("Testing ordinates", epsCmp(vOrds(0,1), ords[0]));

  LC_ASSERT("Testing ordinates", epsCmp(vOrds(1,0), ords[0]));
  LC_ASSERT("Testing ordinates", epsCmp(vOrds(1,1), ords[1]));

  LC_ASSERT("Testing ordinates", epsCmp(vOrds(2,0), ords[1]));
  LC_ASSERT("Testing ordinates", epsCmp(vOrds(2,1), ords[0]));

  LC_ASSERT("Testing ordinates", epsCmp(vOrds(3,0), ords[1]));
  LC_ASSERT("Testing ordinates", epsCmp(vOrds(3,1), ords[1]));

  LC_ASSERT("Testing weights", epsCmp(vWeights[0], weights[0]*weights[0]));
  LC_ASSERT("Testing weights", epsCmp(vWeights[1], weights[0]*weights[1]));
  LC_ASSERT("Testing weights", epsCmp(vWeights[2], weights[1]*weights[0]));
  LC_ASSERT("Testing weights", epsCmp(vWeights[3], weights[1]*weights[1]));

  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(0,0), .6220084679281462, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(0,1), .1666666666666666, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(0,2), .1666666666666666, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(0,3), .04465819873852044, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(1,0), .1666666666666666, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(1,1), .6220084679281462, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(1,2), .04465819873852044, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(1,3), .1666666666666666, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(2,0), .1666666666666666, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(2,1), .04465819873852044, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(2,2), .6220084679281462, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(2,3), .1666666666666666, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(3,0), .04465819873852044, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(3,1), .1666666666666666, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(3,2), .1666666666666666, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(3,3), .6220084679281462, 2e-15));

  unsigned numLowerFaceNodes = basis.getNumSurfLowerNodes(0);
// check quadrature data
  Lucee::Matrix<double> lfInterpMat(numLowerFaceNodes,numTotalNodes), lfOrds(numLowerFaceNodes,2);
  std::vector<double> lfWeights(numLowerFaceNodes);
  basis.getSurfLowerGaussQuadData(0, lfInterpMat, lfOrds, lfWeights);

  LC_ASSERT("Testing surface quadrature ordinates", epsCmp(lfOrds(0,0), -1.0));
  LC_ASSERT("Testing surface quadrature ordinates", epsCmp(lfOrds(0,1), ords[0]));

  LC_ASSERT("Testing surface quadrature ordinates", epsCmp(lfOrds(1,0), -1.0));
  LC_ASSERT("Testing surface quadrature ordinates", epsCmp(lfOrds(1,1), ords[1]));

  LC_ASSERT("Testing surface quadrature weights", epsCmp(lfWeights[0], weights[0]));
  LC_ASSERT("Testing surface quadrature weights", epsCmp(lfWeights[1], weights[1]));

  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(0,0), .7886751345948129, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(0,1), .2113248654051871, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(0,2), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(0,3), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(1,0), .2113248654051871, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(1,1), .7886751345948129, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(1,2), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(1,3), 0.0, 2e-15));

  numLowerFaceNodes = basis.getNumSurfLowerNodes(1);
// check quadrature data
  basis.getSurfLowerGaussQuadData(1, lfInterpMat, lfOrds, lfWeights);

  LC_ASSERT("Testing surface quadrature ordinates", epsCmp(lfOrds(0,0), ords[0]));
  LC_ASSERT("Testing surface quadrature ordinates", epsCmp(lfOrds(0,1), -1.0));

  LC_ASSERT("Testing surface quadrature ordinates", epsCmp(lfOrds(1,0), ords[1]));
  LC_ASSERT("Testing surface quadrature ordinates", epsCmp(lfOrds(1,1), -1.0));

  LC_ASSERT("Testing surface quadrature weights", epsCmp(lfWeights[0], weights[0]));
  LC_ASSERT("Testing surface quadrature weights", epsCmp(lfWeights[1], weights[1]));

  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(0,0), .7886751345948129, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(0,1), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(0,2), .2113248654051871, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(0,3), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(1,0), .2113248654051871, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(1,1), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(1,2), .7886751345948129, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(1,3), 0.0, 2e-15));

  unsigned numUpperFaceNodes = basis.getNumSurfUpperNodes(0);
// check quadrature data
  Lucee::Matrix<double> ufInterpMat(numUpperFaceNodes,numTotalNodes), ufOrds(numUpperFaceNodes,2);
  std::vector<double> ufWeights(numUpperFaceNodes);
  basis.getSurfUpperGaussQuadData(0, ufInterpMat, ufOrds, ufWeights);

  LC_ASSERT("Testing surface quadrature ordinates", epsCmp(ufOrds(0,0), 1.0));
  LC_ASSERT("Testing surface quadrature ordinates", epsCmp(ufOrds(0,1), ords[0]));

  LC_ASSERT("Testing surface quadrature ordinates", epsCmp(ufOrds(1,0), 1.0));
  LC_ASSERT("Testing surface quadrature ordinates", epsCmp(ufOrds(1,1), ords[1]));

  LC_ASSERT("Testing surface quadrature weights", epsCmp(ufWeights[0], weights[0]));
  LC_ASSERT("Testing surface quadrature weights", epsCmp(ufWeights[1], weights[1]));

  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(0,0), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(0,1), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(0,2), .7886751345948129, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(0,3), .2113248654051871, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(1,0), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(1,1), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(1,2), .2113248654051871, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(1,3), .7886751345948129, 2e-15));

  numUpperFaceNodes = basis.getNumSurfUpperNodes(1);
// check quadrature data
  basis.getSurfUpperGaussQuadData(1, ufInterpMat, ufOrds, ufWeights);

  LC_ASSERT("Testing surface quadrature ordinates", epsCmp(ufOrds(0,0), ords[0]));
  LC_ASSERT("Testing surface quadrature ordinates", epsCmp(ufOrds(0,1), 1.0));

  LC_ASSERT("Testing surface quadrature ordinates", epsCmp(ufOrds(1,0), ords[1]));
  LC_ASSERT("Testing surface quadrature ordinates", epsCmp(ufOrds(1,1), 1.0));

  LC_ASSERT("Testing surface quadrature weights", epsCmp(ufWeights[0], weights[0]));
  LC_ASSERT("Testing surface quadrature weights", epsCmp(ufWeights[1], weights[1]));

  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(0,0), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(0,1), .7886751345948129, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(0,2), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(0,3), .2113248654051871, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(1,0), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(1,1), .2113248654051871, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(1,2), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(1,3), .7886751345948129, 2e-15));

// testing nodal weights
  std::vector<double> nodalWeights(basis.getNumNodes());
  basis.getWeights(nodalWeights);

  LC_ASSERT("Testing nodal weights", epsCmp(nodalWeights[0], 1.0));
  LC_ASSERT("Testing nodal weights", epsCmp(nodalWeights[1], 1.0));
  LC_ASSERT("Testing nodal weights", epsCmp(nodalWeights[2], 1.0));
  LC_ASSERT("Testing nodal weights", epsCmp(nodalWeights[3], 1.0));
}

void
test_5_3x3()
{ // 2D, 3x3 nodes
  Lucee::LagrangeTensorBasisCalc<2> basis;
  unsigned nn = 3;
  unsigned numNodes[2];
  numNodes[0] = nn; numNodes[1] = nn;
  basis.calc(Lucee::UNIFORM, numNodes);

  LC_ASSERT("Checking number of nodes", basis.getNumNodes() == nn*nn);
  
// stiffness matrices
  double dx[4] = {2.0, 2.0, 2.0, 2.0};
  Lucee::Matrix<double> stiff(basis.getNumNodes(), basis.getNumNodes());
  basis.getStiffnessMatrix(dx, stiff);

  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(0,0), 28.0/45.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(0,6), (-1.0)/30.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(0,8), (-1.0)/45.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(0,2), (-1.0)/30.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(0,3), (-1.0)/5.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(0,7), 1.0/9.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(0,5), 1.0/9.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(0,1), (-1.0)/5.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(0,4), (-16.0)/45.0, 5e-15));

  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(6,0), (-1.0)/30.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(6,6), 28.0/45.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(6,8), (-1.0)/30.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(6,2), (-1.0)/45.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(6,3), (-1.0)/5.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(6,7), (-1.0)/5.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(6,5), 1.0/9.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(6,1), 1.0/9.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(6,4), (-16.0)/45.0, 5e-15));

  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(8,0), (-1.0)/45.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(8,6), (-1.0)/30.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(8,8), 28.0/45.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(8,2), (-1.0)/30.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(8,3), 1.0/9.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(8,7), (-1.0)/5.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(8,5), (-1.0)/5.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(8,1), 1.0/9.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(8,4), (-16.0)/45.0, 5e-15));

  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(2,0), (-1.0)/30.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(2,6), (-1.0)/45.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(2,8), (-1.0)/30.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(2,2), 28.0/45.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(2,3), 1.0/9.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(2,7), 1.0/9.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(2,5), (-1.0)/5.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(2,1), (-1.0)/5.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(2,4), (-16.0)/45.0, 5e-15));

  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(3,0), (-1.0)/5.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(3,6), (-1.0)/5.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(3,8), 1.0/9.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(3,2), 1.0/9.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(3,3), 88.0/45.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(3,7), (-16.0)/45.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(3,5), 0.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(3,1), (-16.0)/45.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(3,4), (-16.0)/15.0, 5e-15));

  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(7,0), 1.0/9.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(7,6), (-1.0)/5.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(7,8), (-1.0)/5.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(7,2), 1.0/9.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(7,3), (-16.0)/45.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(7,7), 88.0/45.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(7,5), (-16.0)/45.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(7,1), 0.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(7,4), (-16.0)/15.0, 5e-15));

  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(5,0), 1.0/9.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(5,6), 1.0/9.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(5,8), (-1.0)/5.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(5,2), (-1.0)/5.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(5,3), 0.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(5,7), (-16.0)/45.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(5,5), 88.0/45.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(5,1), (-16.0)/45.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(5,4), (-16.0)/15.0, 5e-15));

  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(1,0), (-1.0)/5.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(1,6), 1.0/9.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(1,8), 1.0/9.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(1,2), (-1.0)/5.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(1,3), (-16.0)/45.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(1,7), 0.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(1,5), (-16.0)/45.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(1,1), 88.0/45.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(1,4), (-16.0)/15.0, 5e-15));

  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(4,0), (-16.0)/45.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(4,6), (-16.0)/45.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(4,8), (-16.0)/45.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(4,2), (-16.0)/45.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(4,3), (-16.0)/15.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(4,7), (-16.0)/15.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(4,5), (-16.0)/15.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(4,1), (-16.0)/15.0, 5e-15));
  LC_ASSERT("Testing stiffness matrix", diffCmp(stiff(4,4), 256.0/45.0, 5e-15));

// list of indices
  std::vector<std::vector<int> > eni;
  basis.getExclusiveNodeNdimIndices(eni);
  
  LC_ASSERT("Testing number of exclusive nodes", eni.size() == 4);
  LC_ASSERT("Testing number of exclusive nodes", eni[0].size() == 2);

  LC_ASSERT("Testing exclusive indices", eni[0][0] == 0);
  LC_ASSERT("Testing exclusive indices", eni[0][1] == 0);

  LC_ASSERT("Testing exclusive indices", eni[1][0] == 0);
  LC_ASSERT("Testing exclusive indices", eni[1][1] == 1);

  LC_ASSERT("Testing exclusive indices", eni[2][0] == 1);
  LC_ASSERT("Testing exclusive indices", eni[2][1] == 0);

  LC_ASSERT("Testing exclusive indices", eni[3][0] == 1);
  LC_ASSERT("Testing exclusive indices", eni[3][1] == 1);

  unsigned numTotalNodes = basis.getNumNodes();
// check quadrature data
  Lucee::Matrix<double> vInterpMat(numTotalNodes,numTotalNodes), vOrds(numTotalNodes,2);
  std::vector<double> vWeights(numTotalNodes);
  basis.getGaussQuadData(vInterpMat, vOrds, vWeights);

  double ords[3] = {-std::sqrt(3)/std::sqrt(5), 0.0, std::sqrt(3)/std::sqrt(5)};
  double weights[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};

  unsigned nbasis = 0;
  for (unsigned i=0; i<3; ++i)
    for (unsigned j=0; j<3; ++j)
    {
      LC_ASSERT("Testing ordinates", epsCmp(vOrds(nbasis,0), ords[i]));
      LC_ASSERT("Testing ordinates", epsCmp(vOrds(nbasis,1), ords[j]));
      nbasis++;
    }

  nbasis = 0;
  for (unsigned i=0; i<3; ++i)
    for (unsigned j=0; j<3; ++j)
    {
      LC_ASSERT("Testing weights", epsCmp(vWeights[nbasis], weights[i]*weights[j]));
      nbasis++;
    }

  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(0,0), .4723790007724452, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(0,1), .2749193338482965, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(0,2), -.05999999999999993, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(0,3), .2749193338482966, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(0,4), .1599999999999998, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(0,5), -.03491933384829663, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(0,6), -.05999999999999989, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(0,7), -0.0349193338482966, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(0,8), .007620999227554933, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(1,0), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(1,1), .6872983346207417, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(1,2), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(1,3), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(1,4), .3999999999999997, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(1,5), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(1,6), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(1,7), -.08729833462074155, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(1,8), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(2,0), -.06000000000000006, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(2,1), .2749193338482968, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(2,2), 0.472379000772445, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(2,3), -.03491933384829671, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(2,4), .1599999999999999, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(2,5), .2749193338482965, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(2,6), .007620999227554954, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(2,7), -.03491933384829663, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(2,8), -.05999999999999982, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(3,0), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(3,1), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(3,2), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(3,3), .6872983346207419, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(3,4), .3999999999999997, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(3,5), -.08729833462074163, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(3,6), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(3,7), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(3,8), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(4,0), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(4,1), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(4,2), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(4,3), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(4,4), 1.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(4,5), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(4,6), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(4,7), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(4,8), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(5,0), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(5,1), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(5,2), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(5,3), -.08729833462074182, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(5,4), .4000000000000001, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(5,5), .6872983346207416, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(5,6), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(5,7), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(5,8), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(6,0), -.06000000000000003, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(6,1), -.03491933384829667, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(6,2), .007620999227554949, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(6,3), .2749193338482968, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(6,4), .1599999999999999, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(6,5), -.03491933384829667, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(6,6), 0.472379000772445, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(6,7), .2749193338482965, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(6,8), -.05999999999999991, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(7,0), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(7,1), -.08729833462074174, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(7,2), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(7,3), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(7,4), .4000000000000001, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(7,5), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(7,6), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(7,7), .6872983346207415, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(7,8), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(8,0), .007620999227554969, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(8,1), -.03491933384829671, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(8,2), -.05999999999999994, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(8,3), -.03491933384829674, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(8,4), .1600000000000001, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(8,5), .2749193338482968, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(8,6), -.06000000000000005, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(8,7), .2749193338482967, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(vInterpMat(8,8), .4723790007724449, 2e-15));

  unsigned numLowerFaceNodes = basis.getNumSurfLowerNodes(0);
// check quadrature data
  Lucee::Matrix<double> lfInterpMat(numLowerFaceNodes,numTotalNodes), lfOrds(numLowerFaceNodes,2);
  std::vector<double> lfWeights(numLowerFaceNodes);
  basis.getSurfLowerGaussQuadData(0, lfInterpMat, lfOrds, lfWeights);

  for (unsigned i=0; i<numLowerFaceNodes; ++i)
  {
    LC_ASSERT("Testing surface quadrature ordinates", epsCmp(lfOrds(i,0), -1.0));
    LC_ASSERT("Testing surface quadrature ordinates", epsCmp(lfOrds(i,1), ords[i]));
  }

  for (unsigned i=0; i<numLowerFaceNodes; ++i)
    LC_ASSERT("Testing surface quadrature weights", epsCmp(lfWeights[i], weights[i]));

  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(0,0), .6872983346207417, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(0,1), .3999999999999997, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(0,2), -.08729833462074155, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(0,3), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(0,4), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(0,5), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(0,6), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(0,7), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(0,8), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(1,0), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(1,1), 1.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(1,2), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(1,3), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(1,4), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(1,5), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(1,6), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(1,7), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(1,8), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(2,0), -.08729833462074174, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(2,1), .4000000000000001, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(2,2), .6872983346207415, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(2,3), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(2,4), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(2,5), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(2,6), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(2,7), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(2,8), 0.0, 2e-15));

  basis.getSurfLowerGaussQuadData(1, lfInterpMat, lfOrds, lfWeights);

  for (unsigned i=0; i<numLowerFaceNodes; ++i)
  {
    LC_ASSERT("Testing surface quadrature ordinates", epsCmp(lfOrds(i,0), ords[i]));
    LC_ASSERT("Testing surface quadrature ordinates", epsCmp(lfOrds(i,1), -1.0));
  }

  for (unsigned i=0; i<numLowerFaceNodes; ++i)
    LC_ASSERT("Testing surface quadrature weights", epsCmp(lfWeights[i], weights[i]));

  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(0,0), .6872983346207417, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(0,1), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(0,2), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(0,3), .3999999999999997, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(0,4), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(0,5), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(0,6), -.08729833462074155, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(0,7), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(0,8), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(1,0), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(1,1), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(1,2), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(1,3), 1.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(1,4), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(1,5), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(1,6), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(1,7), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(1,8), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(2,0), -.08729833462074174, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(2,1), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(2,2), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(2,3), .4000000000000001, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(2,4), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(2,5), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(2,6), .6872983346207415, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(2,7), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(lfInterpMat(2,8), 0.0, 2e-15));

  unsigned numUpperFaceNodes = basis.getNumSurfUpperNodes(0);
// check quadrature data
  Lucee::Matrix<double> ufInterpMat(numUpperFaceNodes,numTotalNodes), ufOrds(numUpperFaceNodes,2);
  std::vector<double> ufWeights(numUpperFaceNodes);
  basis.getSurfUpperGaussQuadData(0, ufInterpMat, ufOrds, ufWeights);

  for (unsigned i=0; i<numUpperFaceNodes; ++i)
  {
    LC_ASSERT("Testing surface quadrature ordinates", epsCmp(ufOrds(i,0), 1.0));
    LC_ASSERT("Testing surface quadrature ordinates", epsCmp(ufOrds(i,1), ords[i]));
  }

  for (unsigned i=0; i<numUpperFaceNodes; ++i)
    LC_ASSERT("Testing surface quadrature weights", epsCmp(ufWeights[i], weights[i]));

  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(0,0), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(0,1), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(0,2), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(0,3), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(0,4), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(0,5), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(0,6), .6872983346207417, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(0,7), .3999999999999997, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(0,8), -.08729833462074155, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(1,0), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(1,1), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(1,2), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(1,3), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(1,4), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(1,5), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(1,6), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(1,7), 1.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(1,8), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(2,0), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(2,1), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(2,2), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(2,3), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(2,4), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(2,5), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(2,6), -.08729833462074174, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(2,7), .4000000000000001, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(2,8), .6872983346207415, 2e-15));

  basis.getSurfUpperGaussQuadData(1, ufInterpMat, ufOrds, ufWeights);

  for (unsigned i=0; i<numUpperFaceNodes; ++i)
  {
    LC_ASSERT("Testing surface quadrature ordinates", epsCmp(ufOrds(i,0), ords[i]));
    LC_ASSERT("Testing surface quadrature ordinates", epsCmp(ufOrds(i,1), 1.0));
  }

  for (unsigned i=0; i<numUpperFaceNodes; ++i)
    LC_ASSERT("Testing surface quadrature weights", epsCmp(ufWeights[i], weights[i]));

// testing nodal weights
  std::vector<double> nodalWeights(basis.getNumNodes());
  basis.getWeights(nodalWeights);

  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(0,0), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(0,1), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(0,2), .6872983346207417, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(0,3), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(0,4), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(0,5), .3999999999999997, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(0,6), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(0,7), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(0,8), -.08729833462074155, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(1,0), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(1,1), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(1,2), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(1,3), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(1,4), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(1,5), 1.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(1,6), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(1,7), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(1,8), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(2,0), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(2,1), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(2,2), -.08729833462074174, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(2,3), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(2,4), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(2,5), .4000000000000001, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(2,6), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(2,7), 0.0, 2e-15));
  LC_ASSERT("Testing interpolation matrix", diffCmp(ufInterpMat(2,8), .6872983346207415, 2e-15));

  LC_ASSERT("Testing nodal weights", diffCmp(nodalWeights[0], 1.0/9.0, 2e-15));
  LC_ASSERT("Testing nodal weights", diffCmp(nodalWeights[1], 4.0/9.0, 2e-15));
  LC_ASSERT("Testing nodal weights", diffCmp(nodalWeights[2], 1.0/9.0, 2e-15));
  LC_ASSERT("Testing nodal weights", diffCmp(nodalWeights[3], 4.0/9.0, 2e-15));
  LC_ASSERT("Testing nodal weights", diffCmp(nodalWeights[4], 16.0/9.0, 2e-15));
  LC_ASSERT("Testing nodal weights", diffCmp(nodalWeights[5], 4.0/9.0, 2e-15));
  LC_ASSERT("Testing nodal weights", diffCmp(nodalWeights[6], 1.0/9.0, 2e-15));
  LC_ASSERT("Testing nodal weights", diffCmp(nodalWeights[7], 4.0/9.0, 2e-15));
  LC_ASSERT("Testing nodal weights", diffCmp(nodalWeights[8], 1.0/9.0, 2e-15));
}

void
test_6()
{ // 2D, 4x4 nodes
  Lucee::LagrangeTensorBasisCalc<2> basis;
  unsigned nn = 4;
  unsigned numNodes[2];
  numNodes[0] = nn; numNodes[1] = nn;
  basis.calc(Lucee::GAUSSIAN, numNodes);

  LC_ASSERT("Checking number of nodes", basis.getNumNodes() == nn*nn);
  
  double x1 = std::sqrt((3-2*sqrt(6./5.))/7);
  double x2 = std::sqrt((3+2*sqrt(6./5.))/7);

  double xn[2];
  double loc1[4] = {-x2, -x1, x1, x2};

  unsigned nIdx = 0;
  for (unsigned i=0; i<nn; ++i)
  {
    for (unsigned j=0; j<nn; ++j)
    {
      basis.fillWithNodeCoordinate(nIdx, xn);

      LC_ASSERT("Testing node coordinate", epsCmp(xn[0], loc1[i]));
      LC_ASSERT("Testing node coordinate", epsCmp(xn[1], loc1[j]));

      nIdx += 1;
    }
  }

  nIdx = 0;
  for (unsigned i=0; i<nn; ++i)
  {
    for (unsigned j=0; j<nn; ++j)
    {
      basis.fillWithNodeCoordinate(nIdx, xn);

      for (unsigned b=0; b<nn*nn; ++b)
      {
        if (nIdx == b)
          LC_ASSERT("Testing basis function eval", epsCmp(basis.evalBasis(b, xn), 1.0));
        else
          LC_ASSERT("Testing basis function eval", epsCmp(basis.evalBasis(b, xn), 0.0));
      }

      nIdx += 1;
    }
  }

}

void
test_6_u()
{ // 2D, 4x4 nodes
  Lucee::LagrangeTensorBasisCalc<2> basis;
  unsigned nn = 4;
  unsigned numNodes[2];
  numNodes[0] = nn; numNodes[1] = nn;
  basis.calc(Lucee::UNIFORM, numNodes);

  LC_ASSERT("Checking number of nodes", basis.getNumNodes() == nn*nn);
  
  double x1 = 1-2.0/3.0;
  double x2 = 1.0;

  double xn[2];
  double loc1[4] = {-x2, -x1, x1, x2};

  unsigned nIdx = 0;
  for (unsigned i=0; i<nn; ++i)
  {
    for (unsigned j=0; j<nn; ++j)
    {
      basis.fillWithNodeCoordinate(nIdx, xn);

      LC_ASSERT("Testing node coordinate", epsCmp(xn[0], loc1[i]));
      LC_ASSERT("Testing node coordinate", epsCmp(xn[1], loc1[j]));

      nIdx += 1;
    }
  }

  nIdx = 0;
  for (unsigned i=0; i<nn; ++i)
  {
    for (unsigned j=0; j<nn; ++j)
    {
      basis.fillWithNodeCoordinate(nIdx, xn);

      for (unsigned b=0; b<nn*nn; ++b)
      {
        if (nIdx == b)
          LC_ASSERT("Testing basis function eval", epsCmp(basis.evalBasis(b, xn), 1.0));
        else
          LC_ASSERT("Testing basis function eval", epsCmp(basis.evalBasis(b, xn), 0.0));
      }

      nIdx += 1;
    }
  }

// node lists
  std::vector<int> en;
  basis.getExclusiveNodeIndices(en);

  LC_ASSERT("Testing number of exclusive nodes", en.size() == 9);
  LC_ASSERT("Testing exclusive indices", en[0] == 0);
  LC_ASSERT("Testing exclusive indices", en[1] == 1);
  LC_ASSERT("Testing exclusive indices", en[2] == 2);
  LC_ASSERT("Testing exclusive indices", en[3] == 4);
  LC_ASSERT("Testing exclusive indices", en[4] == 5);
  LC_ASSERT("Testing exclusive indices", en[5] == 6);
  LC_ASSERT("Testing exclusive indices", en[6] == 8);
  LC_ASSERT("Testing exclusive indices", en[7] == 9);
  LC_ASSERT("Testing exclusive indices", en[8] == 10);

  std::vector<int> fn;

  LC_ASSERT("Testing number of lower surface nodes", basis.getNumSurfLowerNodes(0) == 4);
  basis.getSurfLowerNodeNums(0, fn);
  LC_ASSERT("Tesing lower node numbers", fn[0] == 0);
  LC_ASSERT("Tesing lower node numbers", fn[1] == 1);
  LC_ASSERT("Tesing lower node numbers", fn[2] == 2);
  LC_ASSERT("Tesing lower node numbers", fn[3] == 3);

  LC_ASSERT("Testing number of lower surface nodes", basis.getNumSurfLowerNodes(1) == 4);
  basis.getSurfLowerNodeNums(1, fn);
  LC_ASSERT("Tesing lower node numbers", fn[0] == 0);
  LC_ASSERT("Tesing lower node numbers", fn[1] == 4);
  LC_ASSERT("Tesing lower node numbers", fn[2] == 8);
  LC_ASSERT("Tesing lower node numbers", fn[3] == 12);

  LC_ASSERT("Testing number of upper surface nodes", basis.getNumSurfUpperNodes(0) == 4);
  basis.getSurfUpperNodeNums(0, fn);
  LC_ASSERT("Tesing lower node numbers", fn[0] == 12);
  LC_ASSERT("Tesing lower node numbers", fn[1] == 13);
  LC_ASSERT("Tesing lower node numbers", fn[2] == 14);
  LC_ASSERT("Tesing lower node numbers", fn[3] == 15);

  LC_ASSERT("Testing number of upper surface nodes", basis.getNumSurfUpperNodes(1) == 4);
  basis.getSurfUpperNodeNums(1, fn);
  LC_ASSERT("Tesing lower node numbers", fn[0] == 3);
  LC_ASSERT("Tesing lower node numbers", fn[1] == 7);
  LC_ASSERT("Tesing lower node numbers", fn[2] == 11);
  LC_ASSERT("Tesing lower node numbers", fn[3] == 15);

// list of indices
  std::vector<std::vector<int> > eni;
  basis.getExclusiveNodeNdimIndices(eni);

  LC_ASSERT("Testing number of exclusive nodes", eni.size() == 9);
  LC_ASSERT("Testing number of exclusive nodes", eni[0].size() == 2);

  LC_ASSERT("Testing exclusive indices", eni[0][0] == 0);
  LC_ASSERT("Testing exclusive indices", eni[0][1] == 0);

  LC_ASSERT("Testing exclusive indices", eni[1][0] == 0);
  LC_ASSERT("Testing exclusive indices", eni[1][1] == 1);

  LC_ASSERT("Testing exclusive indices", eni[2][0] == 0);
  LC_ASSERT("Testing exclusive indices", eni[2][1] == 2);

  LC_ASSERT("Testing exclusive indices", eni[3][0] == 1);
  LC_ASSERT("Testing exclusive indices", eni[3][1] == 0);

  LC_ASSERT("Testing exclusive indices", eni[4][0] == 1);
  LC_ASSERT("Testing exclusive indices", eni[4][1] == 1);

  LC_ASSERT("Testing exclusive indices", eni[5][0] == 1);
  LC_ASSERT("Testing exclusive indices", eni[5][1] == 2);

  LC_ASSERT("Testing exclusive indices", eni[6][0] == 2);
  LC_ASSERT("Testing exclusive indices", eni[6][1] == 0);

  LC_ASSERT("Testing exclusive indices", eni[7][0] == 2);
  LC_ASSERT("Testing exclusive indices", eni[7][1] == 1);

  LC_ASSERT("Testing exclusive indices", eni[8][0] == 2);
  LC_ASSERT("Testing exclusive indices", eni[8][1] == 2);
}

void
test_7_2()
{  // 3D, 2x2x2 nodes
  Lucee::LagrangeTensorBasisCalc<3> basis;
  unsigned nn = 2;
  unsigned numNodes[3];
  numNodes[0] = nn; numNodes[1] = nn; numNodes[2] = nn;
  basis.calc(Lucee::UNIFORM, numNodes);

  LC_ASSERT("Checking number of nodes", basis.getNumNodes() == nn*nn*nn);
  
  double xn[3];
  double loc1[4] = {-1.0, 1.0};

  unsigned nIdx = 0;
  for (unsigned i=0; i<nn; ++i)
  {
    for (unsigned j=0; j<nn; ++j)
    {
      for (unsigned k=0; k<nn; ++k)
      {
        basis.fillWithNodeCoordinate(nIdx, xn);
        
        LC_ASSERT("Testing node coordinate", epsCmp(xn[0], loc1[i]));
        LC_ASSERT("Testing node coordinate", epsCmp(xn[1], loc1[j]));
        LC_ASSERT("Testing node coordinate", epsCmp(xn[2], loc1[k]));
        
        nIdx += 1;
      }
    }
  }

  nIdx = 0;
  for (unsigned i=0; i<nn; ++i)
  {
    for (unsigned j=0; j<nn; ++j)
    {
      for (unsigned k=0; k<nn; ++k)
      {
        basis.fillWithNodeCoordinate(nIdx, xn);

        for (unsigned b=0; b<nn*nn*nn; ++b)
        {
          if (nIdx == b)
            LC_ASSERT("Testing basis function eval", epsCmp(basis.evalBasis(b, xn), 1.0));
          else
            LC_ASSERT("Testing basis function eval", epsCmp(basis.evalBasis(b, xn), 0.0));
        }

        nIdx += 1;
      }
    }
  }

// node lists
  std::vector<int> en;
  basis.getExclusiveNodeIndices(en);

  LC_ASSERT("Testing number of exclusive nodes", en.size() == 1);
  LC_ASSERT("Testing exclusive indices", en[0] == 0);

  std::vector<int> fn;

  LC_ASSERT("Testing number of lower surface nodes", basis.getNumSurfLowerNodes(0) == 4);
  basis.getSurfLowerNodeNums(0, fn);
  LC_ASSERT("Tesing lower node numbers", fn[0] == 0);
  LC_ASSERT("Tesing lower node numbers", fn[1] == 1);
  LC_ASSERT("Tesing lower node numbers", fn[2] == 2);
  LC_ASSERT("Tesing lower node numbers", fn[3] == 3);

  LC_ASSERT("Testing number of lower surface nodes", basis.getNumSurfLowerNodes(1) == 4);
  basis.getSurfLowerNodeNums(1, fn);
  LC_ASSERT("Tesing lower node numbers", fn[0] == 0);
  LC_ASSERT("Tesing lower node numbers", fn[1] == 1);
  LC_ASSERT("Tesing lower node numbers", fn[2] == 4);
  LC_ASSERT("Tesing lower node numbers", fn[3] == 5);

  LC_ASSERT("Testing number of lower surface nodes", basis.getNumSurfLowerNodes(2) == 4);
  basis.getSurfLowerNodeNums(2, fn);
  LC_ASSERT("Tesing lower node numbers", fn[0] == 0);
  LC_ASSERT("Tesing lower node numbers", fn[1] == 2);
  LC_ASSERT("Tesing lower node numbers", fn[2] == 4);
  LC_ASSERT("Tesing lower node numbers", fn[3] == 6);

  LC_ASSERT("Testing number of upper surface nodes", basis.getNumSurfUpperNodes(0) == 4);
  basis.getSurfUpperNodeNums(0, fn);
  LC_ASSERT("Tesing lower node numbers", fn[0] == 4);
  LC_ASSERT("Tesing lower node numbers", fn[1] == 5);
  LC_ASSERT("Tesing lower node numbers", fn[2] == 6);
  LC_ASSERT("Tesing lower node numbers", fn[3] == 7);

  LC_ASSERT("Testing number of upper surface nodes", basis.getNumSurfUpperNodes(1) == 4);
  basis.getSurfUpperNodeNums(1, fn);
  LC_ASSERT("Tesing lower node numbers", fn[0] == 2);
  LC_ASSERT("Tesing lower node numbers", fn[1] == 3);
  LC_ASSERT("Tesing lower node numbers", fn[2] == 6);
  LC_ASSERT("Tesing lower node numbers", fn[3] == 7);

  LC_ASSERT("Testing number of upper surface nodes", basis.getNumSurfUpperNodes(1) == 4);
  basis.getSurfUpperNodeNums(2, fn);
  LC_ASSERT("Tesing lower node numbers", fn[0] == 1);
  LC_ASSERT("Tesing lower node numbers", fn[1] == 3);
  LC_ASSERT("Tesing lower node numbers", fn[2] == 5);
  LC_ASSERT("Tesing lower node numbers", fn[3] == 7);
}

void
test_7()
{  // 3D, 4x4x4 nodes
  Lucee::LagrangeTensorBasisCalc<3> basis;
  unsigned nn = 4;
  unsigned numNodes[3];
  numNodes[0] = nn; numNodes[1] = nn; numNodes[2] = nn;
  basis.calc(Lucee::GAUSSIAN, numNodes);

  LC_ASSERT("Checking number of nodes", basis.getNumNodes() == nn*nn*nn);
  
  double x1 = std::sqrt((3-2*sqrt(6./5.))/7);
  double x2 = std::sqrt((3+2*sqrt(6./5.))/7);

  double xn[3];
  double loc1[4] = {-x2, -x1, x1, x2};

  unsigned nIdx = 0;
  for (unsigned i=0; i<nn; ++i)
  {
    for (unsigned j=0; j<nn; ++j)
    {
      for (unsigned k=0; k<nn; ++k)
      {
        basis.fillWithNodeCoordinate(nIdx, xn);
        
        LC_ASSERT("Testing node coordinate", epsCmp(xn[0], loc1[i]));
        LC_ASSERT("Testing node coordinate", epsCmp(xn[1], loc1[j]));
        LC_ASSERT("Testing node coordinate", epsCmp(xn[2], loc1[k]));
        
        nIdx += 1;
      }
    }
  }

  nIdx = 0;
  for (unsigned i=0; i<nn; ++i)
  {
    for (unsigned j=0; j<nn; ++j)
    {
      for (unsigned k=0; k<nn; ++k)
      {
        basis.fillWithNodeCoordinate(nIdx, xn);

        for (unsigned b=0; b<nn*nn*nn; ++b)
        {
          if (nIdx == b)
            LC_ASSERT("Testing basis function eval", epsCmp(basis.evalBasis(b, xn), 1.0));
          else
            LC_ASSERT("Testing basis function eval", epsCmp(basis.evalBasis(b, xn), 0.0));
        }

        nIdx += 1;
      }
    }
  }
}

void
test_8()
{  // 3D, 4x3x2 nodes
  Lucee::LagrangeTensorBasisCalc<3> basis;
  unsigned nx = 4, ny = 3, nz = 2;
  unsigned numNodes[3];
  numNodes[0] = nx; numNodes[1] = ny; numNodes[2] = nz;
  basis.calc(Lucee::GAUSSIAN, numNodes);

  LC_ASSERT("Checking number of nodes", basis.getNumNodes() == nx*ny*nz);
  
  double x1 = std::sqrt((3-2*sqrt(6./5.))/7);
  double x2 = std::sqrt((3+2*sqrt(6./5.))/7);

  double xn[3];
  double loc1[4] = {-x2, -x1, x1, x2};
  double loc2[3] = {-sqrt(3.0/5.0), 0, sqrt(3.0/5.0)};
  double loc3[2] = {-1/sqrt(3), 1/sqrt(3)};

  unsigned nIdx = 0;
  for (unsigned i=0; i<nx; ++i)
  {
    for (unsigned j=0; j<ny; ++j)
    {
      for (unsigned k=0; k<nz; ++k)
      {
        basis.fillWithNodeCoordinate(nIdx, xn);
        
        LC_ASSERT("Testing node coordinate", epsCmp(xn[0], loc1[i]));
        LC_ASSERT("Testing node coordinate", epsCmp(xn[1], loc2[j]));
        LC_ASSERT("Testing node coordinate", epsCmp(xn[2], loc3[k]));
        
        nIdx += 1;
      }
    }
  }

  nIdx = 0;
  for (unsigned i=0; i<nx; ++i)
  {
    for (unsigned j=0; j<ny; ++j)
    {
      for (unsigned k=0; k<nz; ++k)
      {
        basis.fillWithNodeCoordinate(nIdx, xn);

        for (unsigned b=0; b<nx*ny*nz; ++b)
        {
          if (nIdx == b)
            LC_ASSERT("Testing basis function eval", epsCmp(basis.evalBasis(b, xn), 1.0));
          else
            LC_ASSERT("Testing basis function eval", epsCmp(basis.evalBasis(b, xn), 0.0));
        }

        nIdx += 1;
      }
    }
  }
}

int
main(int argc, char **argv)
{
  LC_BEGIN_TESTS("lclagrangetensorcalc");
  test_0();
  test_1();
  test_2();
  test_2_lobatto();
  test_3();
  test_4();
  test_4_lobatto();
  test_5();
  test_5_3x3();
  test_6();
  test_6_u();
  test_7_2();
  test_7();
  test_8();
  LC_END_TESTS;
}
