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
  static const unsigned SC_LDG = 0; // local DG
  static const unsigned SC_SLDG = 1; // symmetic local DG
  static const unsigned SC_RDG = 2; // recovery DG

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
    if (scm == "LDG")
      schemeType = SC_LDG;
    else if (scm == "SLDG")
      schemeType = SC_SLDG;
    else if (scm == "RDG")
      schemeType = SC_RDG;
    else
    {
      Lucee::Except lce("DGDiffusionUpdater1D::readInput: unrecongnized scheme '");
      lce << scm << "'."
          << " Must be one of 'LDG', 'SLDG' or 'RDG'";
      throw lce;
    }
  }

  void
  DGDiffusionUpdater1D::initialize()
  {
    Lucee::UpdaterIfc::initialize();

// allocate space
    unsigned nlocal = polyOrder+1;
    iMat = Lucee::Matrix<double>(nlocal, nlocal);
    lowerMat = Lucee::Matrix<double>(nlocal, nlocal);
    upperMat = Lucee::Matrix<double>(nlocal, nlocal);

// set matrices
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
    double cfla = alpha*dt/(dxMin*dxMin);
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
// updater computes inpFld + dt*alpha*d^2/dx^x inpFld
      diffOut.copy(inpFld);
      fact = alpha*dt;
    }
    
    Lucee::ConstFieldPtr<double> inpFldPtr = inpFld.createConstPtr();
    Lucee::FieldPtr<double> diffOutPtr = diffOut.createPtr();

    int idx[1], idxL[1];
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

      idxL[0] = idx[0]+1; // cell attached to upper face
      inpFld.setPtr(inpFldPtr, idxL);
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
}
