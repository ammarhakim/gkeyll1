/**
 * @file	LcDGDiffusionUpdater1D.cpp
 *
 * @brief	Updater to evaluate (hyper)diffusion operators using nodal DG
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLinAlgebra.h>
#include <LcDGDiffusionUpdater1D.h>


namespace Lucee
{
// ID needed by registration system
  const char *DGDiffusionUpdater1D::id = "DGDiffusion1D";

// types of supported schemes
  enum { SC_LDG_L, SC_LDG_R, SC_LDG_S, SC_RDG };


  DGDiffusionUpdater1D::DGDiffusionUpdater1D()
  {
  }

  DGDiffusionUpdater1D::~DGDiffusionUpdater1D()
  {
  }

  void
  DGDiffusionUpdater1D::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

// diffusion coefficient
    alpha = tbl.getNumber("diffusionCoeff");
// CFL number to control time-step
    cfl = tbl.getNumber("cfl"); // CFL number
// should only increments be computed?
    onlyIncrement = false;
    if (tbl.hasBool("onlyIncrement"))
      onlyIncrement = tbl.getBool("onlyIncrement");
// polynomial order
    polyOrder = (int) tbl.getNumber("polyOrder");

// type of scheme
    std::string scm = tbl.getString("scheme");
    if (scm == "LDG-L")
      schemeType = SC_LDG_L;
    else if (scm == "LDG-R")
      schemeType = SC_LDG_R;
    else if (scm == "LDG-S")
      schemeType = SC_LDG_S;
    else if (scm == "RDG")
      schemeType = SC_RDG;
    else
    {
      Lucee::Except lce("DGDiffusionUpdater1D::readInput: unrecongnized scheme '");
      lce << scm << "'."
          << " Must be one of 'LDG-L', 'LDG-R', 'LDG-S' or 'RDG'";
      throw lce;
    }
  }

  void
  DGDiffusionUpdater1D::initialize()
  {
    Lucee::UpdaterIfc::initialize();

    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();
    double dx = grid.getDx(0);

// allocate space
    unsigned nlocal = polyOrder+1;
    iMat = Lucee::Matrix<double>(nlocal, nlocal);
    lowerMat = Lucee::Matrix<double>(nlocal, nlocal);
    upperMat = Lucee::Matrix<double>(nlocal, nlocal);

// set matrices
    switch (schemeType)
    {
      case SC_LDG_L:
          calcLDGLStencil(dx);
          break;
      case SC_LDG_R:
          calcLDGRStencil(dx);
          break;
      case SC_LDG_S:
          calcLDGSStencil(dx);
          break;
      case SC_RDG:
          calcRDGStencil(dx);
          break;
      default:
          // can not happen
          ;
    }
  }

  Lucee::UpdaterStatus
  DGDiffusionUpdater1D::update(double t)
  {
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();

    const Lucee::Field<1, double>& inpFld = this->getInp<Lucee::Field<1, double> >(0);
    Lucee::Field<1, double>& diffOut = this->getOut<Lucee::Field<1, double> >(0);

    double dt = t-this->getCurrTime();
    double dxMin = grid.getDx(0);

// check time-step
    double cflm = 1.1*cfl;
    double cfla = 2.0*alpha*dt/(dxMin*dxMin);
    if (cfla>cflm)
      return Lucee::UpdaterStatus(false, dt*cfl/cfla);

    double fact = 1.0;
    if (onlyIncrement)
    {
// if only increments are requested, the updater computes alpha*d^2/dx^x inpFld
      diffOut = 0.0;
      fact = alpha;
    }
    else
    {
// updater computes inpFld + dt*alpha*d^2/dx^2 inpFld
      diffOut.copy(inpFld);
      fact = alpha*dt;
    }
    
    Lucee::ConstFieldPtr<double> inpFldPtr = inpFld.createConstPtr();
    Lucee::FieldPtr<double> diffOutPtr = diffOut.createPtr();

    int idx[1], idxL[1], idxR[1];
    Lucee::Region<1, int> localRgn = grid.getLocalRegion();

    Lucee::RowMajorSequencer<1> seq(localRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      diffOut.setPtr(diffOutPtr, idx);

// add in contribution to cell from current cell
      inpFld.setPtr(inpFldPtr, idx);
      matVec(fact, iMat, &inpFldPtr[0], 1.0, &diffOutPtr[0]);

// add in contribution from cells attached to lower/upper faces
      idxL[0] = idx[0]-1; // cell attached to lower face
      inpFld.setPtr(inpFldPtr, idxL);
      matVec(fact, lowerMat, &inpFldPtr[0], 1.0, &diffOutPtr[0]);

      idxR[0] = idx[0]+1; // cell attached to upper face
      inpFld.setPtr(inpFldPtr, idxR);
      matVec(fact, upperMat, &inpFldPtr[0], 1.0, &diffOutPtr[0]);
    }

    return Lucee::UpdaterStatus(true, dt*cfl/cfla);
  }

  void
  DGDiffusionUpdater1D::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }

  void 
  DGDiffusionUpdater1D::matVec(double m, const Lucee::Matrix<double>& mat,
    const double* vec, double v, double *out)
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

  void
  DGDiffusionUpdater1D::calcLDGLStencil(double dx)
  {
    double dx2 = dx*dx;

    if (polyOrder == 1)
    {
      lowerMat(0,0) = 4/dx2;
      lowerMat(0,1) = 2/dx2;
      lowerMat(1,0) = -12/dx2;
      lowerMat(1,1) = -6/dx2;

      iMat(0,0) = -8/dx2;
      iMat(0,1) = 2/dx2;
      iMat(1,0) = 6/dx2;
      iMat(1,1) = -24/dx2;

      upperMat(0,0) = 4/dx2;
      upperMat(0,1) = -4/dx2;
      upperMat(1,0) = 6/dx2;
      upperMat(1,1) = -6/dx2;
    }
    else
    {
      throw Lucee::Except("DGDiffusionUpdater1D::calcLDGStencil: Only polyOrder 1 is supported!");
    }
  }

  void
  DGDiffusionUpdater1D::calcLDGRStencil(double dx)
  {
    double dx2 = dx*dx;

    if (polyOrder == 1)
    {
      lowerMat(0,0) = 4/dx2;
      lowerMat(0,1) = 4/dx2;
      lowerMat(1,0) = -6/dx2;
      lowerMat(1,1) = -6/dx2;

      iMat(0,0) = -8/dx2;
      iMat(0,1) = -2/dx2;
      iMat(1,0) = -6/dx2;
      iMat(1,1) = -24/dx2;

      upperMat(0,0) = 4/dx2;
      upperMat(0,1) = -2/dx2;
      upperMat(1,0) = 12/dx2;
      upperMat(1,1) = -6/dx2;
    }
    else
    {
      throw Lucee::Except("DGDiffusionUpdater1D::calcLDGStencil: Only polyOrder 1 is supported!");
    }
  }

  void
  DGDiffusionUpdater1D::calcLDGSStencil(double dx)
  {
    double dx2 = dx*dx;

    if (polyOrder == 1)
    {
      lowerMat(0,0) = 4/dx2;
      lowerMat(0,1) = 3/dx2;
      lowerMat(1,0) = -9/dx2;
      lowerMat(1,1) = -6/dx2;

      iMat(0,0) = -8/dx2;
      iMat(0,1) = 0;
      iMat(1,0) = 0;
      iMat(1,1) = -24/dx2;

      upperMat(0,0) = 4/dx2;
      upperMat(0,1) = -3/dx2;
      upperMat(1,0) = 9/dx2;
      upperMat(1,1) = -6/dx2;
    }
    else
    {
      throw Lucee::Except("DGDiffusionUpdater1D::calcSLDGStencil: Only polyOrder 1 is supported!");
    }
  }

  void
  DGDiffusionUpdater1D::calcRDGStencil(double dx)
  {
    double dx2 = dx*dx;

    if (polyOrder == 1)
    {
      lowerMat(0,0) = 9.0/(4.0*dx2);
      lowerMat(0,1) = 5.0/(4.0*dx2);
      lowerMat(1,0) = (-15.0)/(4.0*dx2);
      lowerMat(1,1) = (-7.0)/(4.0*dx2);

      iMat(0,0) = (-9.0)/(2.0*dx2);
      iMat(0,1) = 0;
      iMat(1,0) = 0;
      iMat(1,1) = (-23.0)/(2.0*dx2);

      upperMat(0,0) = 9.0/(4.0*dx2);
      upperMat(0,1) = (-5.0)/(4.0*dx2);
      upperMat(1,0) = 15.0/(4.0*dx2);
      upperMat(1,1) = (-7.0)/(4.0*dx2);
    }
    else
    {
      throw Lucee::Except("DGDiffusionUpdater1D::calcRDGStencil: Only polyOrder 1 is supported!");
    }
  }
}
