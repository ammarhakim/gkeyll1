/**
 * @file	LcSOLBGKCollisionUpdater.cpp
 *
 * @brief	Updater to compute the drag term in the L-B collision operator.
 * df/dt = alpha*d[(v-uIn)f]/dv where uIn is a function of (x,y,z) (drift velocity)
 * Used for 3D2V SOL problem.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLBGKCollisionUpdater.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
// set id for module system
  const char *SOLBGKCollisionUpdater::id = "SOLBGKCollisionUpdater";

  SOLBGKCollisionUpdater::SOLBGKCollisionUpdater()
    : UpdaterIfc()
  {
  }

  void 
  SOLBGKCollisionUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except("SOLBGKCollisionUpdater::readInput: Must specify element to use using 'basis5d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("SOLBGKCollisionUpdater::readInput: Must specify element to use using 'basis3d'");
    
    onlyIncrement = false;
    if (tbl.hasBool("onlyIncrement"))
      onlyIncrement = tbl.getBool("onlyIncrement");

    if (tbl.hasFunction("alpha"))
      fnRef = tbl.getFunctionRef("alpha");
    else
      throw Lucee::Except("SOLBGKCollisionUpdater::readInput: Must supply a collision frequency function as alpha.");
  }

  void 
  SOLBGKCollisionUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();
    // local region to update
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();

    Lucee::RowMajorSequencer<5> seq(localRgn);
    seq.step(); // just to get to first index
    int idx[5];
    seq.fillWithIndex(idx);
    nodalBasis5d->setIndex(idx);
    nodalBasis3d->setIndex(idx[0], idx[1], idx[2]);
    
    int nlocal3d = nodalBasis3d->getNumNodes();
    // Compute mom0Vector for cell-average calculation
    Lucee::Matrix<double> massMatrix3dLucee(nlocal3d, nlocal3d);
    nodalBasis3d->getMassMatrix(massMatrix3dLucee);
    Eigen::MatrixXd massMatrix3d(nlocal3d, nlocal3d);
    copyLuceeToEigen(massMatrix3dLucee, massMatrix3d);
    mom0Vector = massMatrix3d.colwise().sum();
  }

  Lucee::UpdaterStatus 
  SOLBGKCollisionUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<5>& grid
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // Distribution function
    const Lucee::Field<5, double>& fIn = this->getInp<Lucee::Field<5, double> >(0);
    // Maxwellian function
    const Lucee::Field<5, double>& fMaxwellian = this->getInp<Lucee::Field<5, double> >(1);
    // Number density at node
    const Lucee::Field<3, double>& numDensityIn = this->getInp<Lucee::Field<3, double> >(2);
    // Temperature in joules
    const Lucee::Field<3, double>& temperatureIn = this->getInp<Lucee::Field<3, double> >(3);
    // Output distribution function
    Lucee::Field<5, double>& fOut = this->getOut<Lucee::Field<5, double> >(0);

    int nlocal5d = nodalBasis5d->getNumNodes();
    int nlocal3d = nodalBasis3d->getNumNodes();

    double dt = t-this->getCurrTime();

    Lucee::Region<5, int> localRgn = grid.getLocalRegion();
    Lucee::Region<5, int> globalRgn = grid.getGlobalRegion();

    double cfla = 0.0; // maximum CFL number

    Lucee::ConstFieldPtr<double> fInPtr = fIn.createConstPtr();
    Lucee::ConstFieldPtr<double> fMaxwellianPtr = fMaxwellian.createConstPtr();
    Lucee::ConstFieldPtr<double> numDensityInPtr = numDensityIn.createConstPtr();
    Lucee::ConstFieldPtr<double> temperatureInPtr = temperatureIn.createConstPtr();
    Lucee::FieldPtr<double> fOutPtr = fOut.createPtr();

    fOut = 0.0; // use fOut to store increment initially
    int idx[5];
    Lucee::RowMajorSequencer<5> seq(localRgn);
    
    // Get alpha. Need to scale by n/T^(3/2) in this class
    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;
    std::vector<double> resultVector(1);
    evaluateFunction(*L, t, resultVector);
    alpha = resultVector[0];

    while (seq.step())
    {
      seq.fillWithIndex(idx);
      fIn.setPtr(fInPtr, idx);
      fMaxwellian.setPtr(fMaxwellianPtr, idx);
      fOut.setPtr(fOutPtr, idx);
      // Compute alpha scale factor n/Te^(3/2)
      temperatureIn.setPtr(temperatureInPtr, idx[0], idx[1], idx[2]);
      numDensityIn.setPtr(numDensityInPtr, idx[0], idx[1], idx[2]);
      Eigen::VectorXd temperatureVec(nlocal3d);
      Eigen::VectorXd numDensityVec(nlocal3d);
      for (int nodeIndex = 0; nodeIndex < nlocal3d; nodeIndex++)
      {
        temperatureVec(nodeIndex) = temperatureInPtr[nodeIndex];
        numDensityVec(nodeIndex) = numDensityInPtr[nodeIndex];
      }
      grid.setIndex(idx);
      double volume = grid.getDx(0)*grid.getDx(1)*grid.getDx(2);
      double averageTemperature = mom0Vector.dot(temperatureVec)/volume;
      double averageN = mom0Vector.dot(numDensityVec)/volume;
      
      if (onlyIncrement == false)
      {
        for (int nodeIndex = 0; nodeIndex < nlocal5d; nodeIndex++)
          fOutPtr[nodeIndex] = fInPtr[nodeIndex] - dt*alpha*averageN/(averageTemperature*sqrt(averageTemperature))*
            (fInPtr[nodeIndex]-fMaxwellianPtr[nodeIndex]);
      }
      else
      {
        for (int nodeIndex = 0; nodeIndex < nlocal5d; nodeIndex++)
          fOutPtr[nodeIndex] = -alpha*averageN/(averageTemperature*sqrt(averageTemperature))*
            (fInPtr[nodeIndex]-fMaxwellianPtr[nodeIndex]);
      }
    }
    return Lucee::UpdaterStatus();
  }

  void
  SOLBGKCollisionUpdater::declareTypes()
  {
    // Input: Old distribution
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Input: Maxwellian distribution
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Input: <1> number density
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: Temperature in joules
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Output: Updated distribution
    this->appendOutVarType(typeid(Lucee::Field<5, double>));
  }

  void
  SOLBGKCollisionUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
    {
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
      {
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
      }
    }
  }

  void
  SOLBGKCollisionUpdater::evaluateFunction(Lucee::LuaState& L, double tm,
    std::vector<double>& res)
  {
    // push function object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, fnRef);
    // push variables on stack
    lua_pushnumber(L, tm);
    // call function
    if (lua_pcall(L, 1, res.size(), 0) != 0)
    {
      Lucee::Except lce("SOLBGKCollisionUpdater::evaluateFunction: ");
      lce << "Problem evaluating function supplied as 'alpha' "
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
        throw Lucee::Except("SOLBGKCollisionUpdater::evaluateFunction: Return value not a number");
      res[res.size()+i] = lua_tonumber(L, i);
    }
    lua_pop(L, 1);
  }
}
