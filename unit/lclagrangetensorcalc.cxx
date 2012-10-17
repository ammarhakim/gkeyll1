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
  test_3();
  test_4();
  test_5();
  test_6();
  test_7();
  test_8();
  LC_END_TESTS;
}
