/**
 * @file	LcSmoothQuadPhiToC1Updater.cpp
 *
 * @brief	Project piecewise quadriatic phi to C1 basis functions
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLinAlgebra.h>
#include <LcSmoothQuadPhiToC1Updater.h>

namespace Lucee
{
  static const char *SmoothQuadPhiToC1Updaterid = "SmoothQuadPhiToC1";

  static
  double avg(double *nv)
  {
    return 0.5*(nv[0]+nv[2]);
  }

  static
  double slp(double *nv)
  {
    return nv[2]-nv[0];
  }

  SmoothQuadPhiToC1Updater::SmoothQuadPhiToC1Updater()
    : UpdaterIfc()
  {
  }

  void
  SmoothQuadPhiToC1Updater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis");
    else
      throw Lucee::Except("SmoothQuadPhiToC1Updater::readInput: Must specify element to use using 'basis'");
  }

  void
  SmoothQuadPhiToC1Updater::initialize()
  {
    Lucee::UpdaterIfc::initialize();

// get hold of grid
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();

// local region to update
    Lucee::Region<1, int> localRgn = grid.getLocalRegion();
    Lucee::RowMajorSequencer<1> seq(localRgn);
    seq.step();
    int idx[1];
    seq.fillWithIndex(idx);

// set index to first location in grid (this is okay as in this
// updater we are assuming grid is uniform)
    nodalBasis->setIndex(idx);

    unsigned nVol = nodalBasis->getNumGaussNodes();
    unsigned nlocal = nodalBasis->getNumNodes();

// get stiffness matrix
    Lucee::Matrix<double> stiffMatrix(nlocal, nlocal);
    nodalBasis->getGradStiffnessMatrix(0, stiffMatrix);

// calculate differentiation matrix
    diffMatrix = Lucee::Matrix<double>(nlocal, nlocal);
    for (unsigned i=0; i<nlocal; ++i)
      for (unsigned j=0; j<nlocal; ++j)
// diff matrices are computed from transposed stiffness matrices
          diffMatrix(i,j) = stiffMatrix(j,i);

// multiply by inverse of mass matrix
    Lucee::Matrix<double> massMatrix(nlocal, nlocal);
    nodalBasis->getMassMatrix(massMatrix);
    Lucee::solve(massMatrix, diffMatrix);
  }

  Lucee::UpdaterStatus
  SmoothQuadPhiToC1Updater::update(double t)
  {
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();
    const Lucee::Field<1, double>& phiC0 = this->getInp<Lucee::Field<1, double> >(0);
    Lucee::Field<1, double>& phiC1 = this->getOut<Lucee::Field<1, double> >(0);

    Lucee::Region<1, int> localRgn = grid.getLocalRegion();
    Lucee::RowMajorSequencer<1> seq(localRgn);

    unsigned nlocal = nodalBasis->getNumNodes();
    std::vector<double> phiK(nlocal);

    std::vector<double> grad(3*(localRgn.getVolume()+2)); // +2 for ghost cells

// compute gradient of potential (loop is over extended region)
    for (unsigned i=localRgn.getLower(0)-1; i<localRgn.getUpper(0)+1; ++i)
    {
      nodalBasis->setIndex(i);
// extract potential at this location
      nodalBasis->extractFromField(phiC0, phiK);
// compute gradient and store in appropriate location
      matVec(1.0, diffMatrix, phiK, 0.0, &grad[3*i]);
    }

    double alpha = 1.0/3.0; // weight factor
    double a2 = 0.5*alpha;

// smooth gradient (this is actually piece-wise linear, even though
// three local nodes are used). Edge values are stored.
    std::vector<double> smoothedGradEdge(localRgn.getVolume()+1); // one more edge than cells
    for (unsigned i=localRgn.getLower(0); i<localRgn.getUpper(0)+1; ++i)
    {
      smoothedGradEdge[i] = 0.5*(
        avg(&grad[3*i]) - a2*slp(&grad[3*i])
        + avg(&grad[3*(i-1)]) + a2*slp(&grad[3*(i-1)]) );
    }    

    double cumTot = smoothedGradEdge[0]; // value at the left-most edge
    Lucee::FieldPtr<double> ptr = phiC1.createPtr();

// integrate to compute a C1 approximation to potential
    double dx = grid.getDx(0);
    for (unsigned i=localRgn.getLower(0); i<localRgn.getUpper(0)+1; ++i)
    {
      phiC1.setPtr(ptr, i);
      ptr[0] = cumTot;
      ptr[1] = ptr[0] + dx/8*(smoothedGradEdge[i+1] + 3*smoothedGradEdge[i]);
      cumTot += ptr[0] + dx/2*(smoothedGradEdge[i+1] + smoothedGradEdge[i]);
    }        

    return Lucee::UpdaterStatus();
  }

  void
  SmoothQuadPhiToC1Updater::declareTypes()
  {
// takes one input (quadriatic, C0 phi)
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
// returns one output, C1 phi
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }

  void 
  SmoothQuadPhiToC1Updater::matVec(double m, const Lucee::Matrix<double>& mat,
    const std::vector<double>& vec, double v, double *out)
  {
    double tv;
    unsigned rows = mat.numRows(), cols = mat.numColumns();
    for (unsigned i=0; i<rows; ++i)
    {
      tv = 0.0;
      for (unsigned j=0; j<cols; ++j)
        tv += mat(i,j)*vec[j];
      out[i] = m*tv + v*out[i];
    }
  }
}
