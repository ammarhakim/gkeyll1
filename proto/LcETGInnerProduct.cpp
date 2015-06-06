/**
 * @file	LcETGInnerProduct.cpp
 *
 * @brief	Updater to compute the free energy for ETG problem.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGlobals.h>
#include <LcETGInnerProduct.h>
#include <LcMathPhysConstants.h>

// loki includes
#include <loki/Singleton.h>

// eigen include to write to file
#include <unsupported/Eigen/SparseExtra>

namespace Lucee
{
  const char *ETGInnerProduct::id = "ETGInnerProduct";

  ETGInnerProduct::ETGInnerProduct()
    : Lucee::UpdaterIfc()
  {
  }

  void
  ETGInnerProduct::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    // get hold of element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<4> >("basis4d"))
      nodalBasis4d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<4> >("basis4d");
    else
      throw Lucee::Except("ETGInnerProduct::readInput: Must specify element to use using 'basis4d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis2d"))
      nodalBasis2d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis2d");
    else
      throw Lucee::Except("ETGInnerProduct::readInput: Must specify element to use using 'basis2d'");

    // in eV
    if (tbl.hasNumber("adiabaticTemp"))
      bgAdiabaticTemp = tbl.getNumber("adiabaticTemp");
    else
      throw Lucee::Except(
        "ETGInnerProduct::readInput: Must specify temperature of adiabatic species using 'adiabaticTemp'");

    if (tbl.hasNumber("adiabaticDensity"))
      bgAdiabaticDensity = tbl.getNumber("adiabaticDensity");
    else
      throw Lucee::Except(
        "ETGInnerProduct::readInput: Must specify density of adiabatic species using 'adiabaticDensity'");

    // For now, assume Z_i*q_i = elementary charge
    if (tbl.hasNumber("eV"))
      eV = tbl.getNumber("eV");
    else
      throw Lucee::Except(
        "ETGInnerProduct::readInput: Must specify conversion of joules to eV using 'eV'");

    if (tbl.hasNumber("kineticMass"))
      kineticMass = tbl.getNumber("kineticMass");
    else
      throw Lucee::Except(
        "ETGInnerProduct::readInput: Must specify mass of kinetic species using 'kineticMass'");

    if (tbl.hasNumber("totalNodes"))
      totalNodes = tbl.getNumber("totalNodes");
    else
      throw Lucee::Except(
        "ETGInnerProduct::readInput: Must specify total number of nodes in system using 'totalNodes'");

    // get function to evaluate
    fnRef = tbl.getFunctionRef("evaluate");
  }

  void
  ETGInnerProduct::initialize()
  {
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<4>& grid 
      = this->getGrid<Lucee::StructuredGridBase<4> >();

    // local region to update
    Lucee::Region<4, int> localRgn = grid.getLocalRegion();

    // set index to first location in grid (this is okay as in this
    // updater we are assuming grid is uniform)
    Lucee::RowMajorSequencer<4> seq(localRgn);
    seq.step();
    int idx[4];
    seq.fillWithIndex(idx);
    nodalBasis4d->setIndex(idx);
    nodalBasis2d->setIndex(idx[0], idx[1]);
    
    int nlocal2d = nodalBasis2d->getNumNodes();
    int nlocal4d = nodalBasis4d->getNumNodes();

    // get data needed for Gaussian quadrature (4D)
    int nVolQuad4d = nodalBasis4d->getNumGaussNodes();
    std::vector<double> volWeights4d(nVolQuad4d);
    Lucee::Matrix<double> tempVolQuad4d(nVolQuad4d, nlocal4d);
    Lucee::Matrix<double> tempVolCoords4d(nVolQuad4d, 4);
    volQuad4d.reset(nVolQuad4d, nlocal4d);

    nodalBasis4d->getGaussQuadData(tempVolQuad4d, tempVolCoords4d, volWeights4d);
    for (int volIndex = 0; volIndex < nVolQuad4d; volIndex++)
      volQuad4d.weights(volIndex) = volWeights4d[volIndex];
    
    copyLuceeToEigen(tempVolQuad4d, volQuad4d.interpMat);

    // get data needed for Gaussian quadrature (2D)
    int nVolQuad2d = nodalBasis2d->getNumGaussNodes();
    std::vector<double> volWeights2d(nVolQuad2d);
    Lucee::Matrix<double> tempVolQuad2d(nVolQuad2d, nlocal2d);
    Lucee::Matrix<double> tempVolCoords2d(nVolQuad2d, 3);
    volQuad2d.reset(nVolQuad2d, nlocal2d);

    nodalBasis2d->getGaussQuadData(tempVolQuad2d, tempVolCoords2d, volWeights2d);
    for (int volIndex = 0; volIndex < nVolQuad2d; volIndex++)
      volQuad2d.weights(volIndex) = volWeights2d[volIndex];
    
    copyLuceeToEigen(tempVolQuad2d, volQuad2d.interpMat);
  }

  Lucee::UpdaterStatus
  ETGInnerProduct::update(double t)
  {
    const Lucee::StructuredGridBase<4>& grid 
      = this->getGrid<Lucee::StructuredGridBase<4> >();

    const Lucee::Field<2, double>& bFieldIn = this->getInp<Lucee::Field<2, double> >(0);
    const Lucee::Field<2, double>& bgKineticTempIn = this->getInp<Lucee::Field<2, double> >(1);
    const Lucee::Field<4, double>& bgDistFIn = this->getInp<Lucee::Field<4, double> >(2);
    const Lucee::Field<2, double>& gPotentialIn = this->getInp<Lucee::Field<2, double> >(3);
    const Lucee::Field<2, double>& hPotentialIn = this->getInp<Lucee::Field<2, double> >(4);
    const Lucee::Field<4, double>& gDistFIn = this->getInp<Lucee::Field<4, double> >(5);
    const Lucee::Field<4, double>& hDistFIn = this->getInp<Lucee::Field<4, double> >(6);

    // iterators into fields
    Lucee::ConstFieldPtr<double> bFieldPtr = bFieldIn.createConstPtr();
    Lucee::ConstFieldPtr<double> bgKineticTempPtr = bgKineticTempIn.createConstPtr();
    Lucee::ConstFieldPtr<double> bgDistFPtr = bgDistFIn.createConstPtr();
    Lucee::ConstFieldPtr<double> gPotentialPtr = gPotentialIn.createConstPtr();
    Lucee::ConstFieldPtr<double> hPotentialPtr = hPotentialIn.createConstPtr();
    Lucee::ConstFieldPtr<double> gDistFPtr = gDistFIn.createConstPtr();
    Lucee::ConstFieldPtr<double> hDistFPtr = hDistFIn.createConstPtr();

    Lucee::Region<4, int> localRgn = grid.getLocalRegion();
    Lucee::Region<4, int> globalRgn = grid.getGlobalRegion();
    Lucee::RowMajorSequencer<4> seq(localRgn);

    unsigned nlocal4d = nodalBasis4d->getNumNodes();
    unsigned nlocal2d = nodalBasis2d->getNumNodes();

    int nVolQuad4d = nodalBasis4d->getNumGaussNodes();
    int nVolQuad2d = nodalBasis2d->getNumGaussNodes();
    int idx[4];

    double LX = (globalRgn.getUpper(0) - globalRgn.getLower(0))*grid.getDx(0);
    double LY = (globalRgn.getUpper(1) - globalRgn.getLower(1))*grid.getDx(1);

    // Get write location
    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;
    std::vector<double> res(2);
    evaluateFunction(*L, t, res);

    double localInt = 0.0;
    // loop, performing integration
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      // set index into element basis
      nodalBasis4d->setIndex(idx);
      nodalBasis2d->setIndex(idx[0], idx[1]);

      bFieldIn.setPtr(bFieldPtr, idx[0], idx[1]);
      bgKineticTempIn.setPtr(bgKineticTempPtr, idx[0], idx[1]);
      bgDistFIn.setPtr(bgDistFPtr, idx);
      gPotentialIn.setPtr(gPotentialPtr, idx[0], idx[1]);
      hPotentialIn.setPtr(hPotentialPtr, idx[0], idx[1]);
      gDistFIn.setPtr(gDistFPtr, idx);
      hDistFIn.setPtr(hDistFPtr, idx);
      
      Eigen::VectorXd bFieldVec(nlocal2d);
      Eigen::VectorXd bgKineticTempVec(nlocal2d);
      Eigen::VectorXd bgDistFVec(nlocal4d);
      Eigen::VectorXd gPotentialVec(nlocal2d);
      Eigen::VectorXd hPotentialVec(nlocal2d);
      Eigen::VectorXd gDistFVec(nlocal4d);
      Eigen::VectorXd hDistFVec(nlocal4d);

      for (int i = 0; i < nlocal2d; i++)
      {
        bFieldVec(i) = bFieldPtr[i];
        bgKineticTempVec(i) = bgKineticTempPtr[i];
        gPotentialVec(i) = gPotentialPtr[i];
        hPotentialVec(i) = hPotentialPtr[i];
      }

      for (int i = 0; i < nlocal4d; i++)
      {
        bgDistFVec(i) = bgDistFPtr[i];
        gDistFVec(i) = gDistFPtr[i];
        hDistFVec(i) = hDistFPtr[i];
      }

      // Interpolate data to quadrature points
      Eigen::VectorXd bFieldAtQuad = volQuad2d.interpMat*bFieldVec;
      Eigen::VectorXd bgKineticTempAtQuad = volQuad2d.interpMat*bgKineticTempVec;
      Eigen::VectorXd bgDistFAtQuad = volQuad4d.interpMat*bgDistFVec;
      Eigen::VectorXd gPotentialAtQuad = volQuad2d.interpMat*gPotentialVec;
      Eigen::VectorXd hPotentialAtQuad = volQuad2d.interpMat*hPotentialVec;
      Eigen::VectorXd gDistFAtQuad = volQuad4d.interpMat*gDistFVec;
      Eigen::VectorXd hDistFAtQuad = volQuad4d.interpMat*hDistFVec;

      // perform quadrature
      for (int i = 0; i < nVolQuad4d; i++)
      {
        localInt += volQuad4d.weights[i]/(LX*LY)*( 2*Lucee::PI*bFieldAtQuad(i % nVolQuad2d)/kineticMass*
          bgKineticTempAtQuad(i % nVolQuad2d)*gDistFAtQuad(i)*hDistFAtQuad(i)/(2*bgDistFAtQuad(i)) +
          bgAdiabaticDensity/(2*bgAdiabaticTemp)*
          gPotentialAtQuad(i % nVolQuad2d)*hPotentialAtQuad(i % nVolQuad2d) );
      }
    }

    double volInt = localInt;
    // get hold of comm pointer to do all parallel messaging
    TxCommBase *comm = this->getComm();
    comm->allreduce(1, &localInt, &volInt, TX_SUM);
    
    int rowIndex = (int) res[0];
    int colIndex = (int) res[1];
 
    if (volInt != 0.0)
    {
      //std::cout << "Nonzero element at (" << rowIndex << "," << colIndex << ")" << std::endl;
      //std::cout << "Value is = " << volInt << std::endl;
      //innerProductMatrix.insert(rowIndex, colIndex) = volInt;
      //if (rowIndex != colIndex)
      //  innerProductMatrix.insert(colIndex, rowIndex) = volInt;
      tripletList.push_back(Eigen::Triplet<double>(rowIndex,colIndex,volInt));
      //if (rowIndex != colIndex)
      //  tripletList.push_back(Eigen::Triplet<double>(colIndex,rowIndex,volInt));
    }

    // check if we need to create and output matrix
    if (rowIndex == totalNodes-1 && colIndex == totalNodes-1)
    {
      Eigen::SparseMatrix<double> innerProductMatrix(totalNodes, totalNodes);
      innerProductMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
      Eigen::saveMarket(innerProductMatrix,"innerProductMatrix.mtx");
      // Perform a Cholesky decomposition of innerProductMatrix
      //Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > choleskyDecomp(innerProductMatrix);
      // Retrieve U (A = U^* U)
      //Eigen::MatrixXd choleskyResult = choleskyDecomp.matrixU();
      //Eigen::saveMarket(choleskyResult,"innerProductMatrixDecomp.mtx");
    }

    return Lucee::UpdaterStatus();
  }

  void
  ETGInnerProduct::declareTypes()
  {
    // Magnetic field B
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    // Background kinetic temperature (in eV)
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    // Background kinetic distribution function
    this->appendInpVarType(typeid(Lucee::Field<4, double>));
    // Smoothed potential from distribution function g
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    // Smoothed potential from distribution function h
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    // Entire distribution function g
    this->appendInpVarType(typeid(Lucee::Field<4, double>));
    // Entire distribution function h
    this->appendInpVarType(typeid(Lucee::Field<4, double>));
  }

  void
  ETGInnerProduct::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }

  void
  ETGInnerProduct::evaluateFunction(Lucee::LuaState& L, double tm,
    std::vector<double>& res)
  {
    // push function object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, fnRef);
    lua_pushnumber(L, tm);
    // call function
    if (lua_pcall(L, 1, res.size(), 0) != 0)
    {
      Lucee::Except lce("ETGInnerProduct::evaluateFunction: ");
      lce << "Problem evaluating function supplied as 'evaluate' "
          << std::endl;
      std::string err(lua_tostring(L, -1));
      lua_pop(L, 1);
      lce << "[" << err << "]";
      throw lce;
    }
    // fetch results
    for (int i=-res.size(); i<0; ++i)
    {
      if (!lua_isnumber(L, i))
        throw Lucee::Except("ETGInnerProduct::evaluateFunction: Return value not a number");
      res[res.size()+i] = lua_tonumber(L, i);
    }
    lua_pop(L, 1);
  }
}
