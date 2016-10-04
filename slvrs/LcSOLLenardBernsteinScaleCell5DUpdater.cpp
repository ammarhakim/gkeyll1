/**
 * @file	LcSOLLenardBernsteinScaleCell5DUpdater.cpp
 *
 * @brief	Accumulates correct amount of an input distribution function (assumed to be a diffusion term)
 * to an existing distribution so that the desired amount of energy in each cell is achieved
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLLenardBernsteinScaleCell5DUpdater.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  const char *SOLLenardBernsteinScaleCell5DUpdater::id = "SOLLenardBernsteinScaleCell5D";

  SOLLenardBernsteinScaleCell5DUpdater::SOLLenardBernsteinScaleCell5DUpdater()
  {
  }

  void
  SOLLenardBernsteinScaleCell5DUpdater::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except("SOLLenardBernsteinScaleCell5DUpdater::readInput: Must specify element to use using 'basis5d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("SOLLenardBernsteinScaleCell5DUpdater::readInput: Must specify element to use using 'basis3d'");

    // CFL number to control time-step
    cfl = tbl.getNumber("cfl"); // CFL number

    if (tbl.hasNumber("speciesMass"))
      speciesMass = tbl.getNumber("speciesMass");
    else
      throw Lucee::Except("SOLLenardBernsteinScaleCell5DUpdater::readInput: Must specify speciesMass");

    if (tbl.hasFunction("alpha"))
      fnRef = tbl.getFunctionRef("alpha");
    else
      throw Lucee::Except("SOLLenardBernsteinScaleCell5DUpdater::readInput: Must supply a collision frequency function as alpha.");
  }

  void
  SOLLenardBernsteinScaleCell5DUpdater::initialize()
  {
    UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus
  SOLLenardBernsteinScaleCell5DUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // Distribution function after diffusion step
    const Lucee::Field<5, double>& distfDiff = this->getInp<Lucee::Field<5, double> >(0);
    // Energy in each cell after drag step
    const Lucee::Field<3, double>& energyDragIn = this->getInp<Lucee::Field<3, double> >(1);
    // Energy in each cell after diffusion step
    const Lucee::Field<3, double>& energyDiffIn = this->getInp<Lucee::Field<3, double> >(2);
    // Temperature in joules
    const Lucee::Field<3, double>& temperatureIn = this->getInp<Lucee::Field<3, double> >(3);
    // Magnetic field
    const Lucee::Field<3, double>& bFieldIn = this->getInp<Lucee::Field<3, double> >(4);
    // Dimensionally correct number density from weighted moment calculation
    const Lucee::Field<3, double>& numDensityIn = this->getInp<Lucee::Field<3, double> >(5);
    // Distribution function with drag already applied and diffusion to be added on in update()
    Lucee::Field<5, double>& distfOut = this->getOut<Lucee::Field<5, double> >(0);
    
    Lucee::ConstFieldPtr<double> distfDiffPtr   = distfDiff.createConstPtr();
    Lucee::ConstFieldPtr<double> energyDragInPtr = energyDragIn.createConstPtr();
    Lucee::ConstFieldPtr<double> energyDiffInPtr = energyDiffIn.createConstPtr();
    Lucee::ConstFieldPtr<double> temperatureInPtr = temperatureIn.createConstPtr();
    Lucee::ConstFieldPtr<double> bFieldInPtr = bFieldIn.createConstPtr();
    Lucee::ConstFieldPtr<double> numDensityInPtr = numDensityIn.createConstPtr();

    // Remember not to clear distfOut
    Lucee::FieldPtr<double> distfOutPtr = distfOut.createPtr(); // Output pointer

    unsigned nlocal5d = nodalBasis5d->getNumNodes();
    unsigned nlocal3d = nodalBasis3d->getNumNodes();
    double dt = t-this->getCurrTime();

    double cellCentroid[5];
    int idx[5];

    Lucee::Region<5, int> globalRgn = grid.getGlobalRegion();
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();

    // check time-step
    double cflm = 1.1*cfl;
    double cfla = 0.0;
    // Get alpha. Need to scale by n/T^(3/2) in this class
    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;
    std::vector<double> resultVector(1);
    evaluateFunction(*L, t, resultVector);
    double alpha = resultVector[0];

    Lucee::RowMajorSequencer<5> seq(localRgn);

    while (seq.step())
    {
      seq.fillWithIndex(idx);
      grid.setIndex(idx);
      grid.getCentroid(cellCentroid);

      energyDragIn.setPtr(energyDragInPtr, idx[0], idx[1], idx[2]);
      energyDiffIn.setPtr(energyDiffInPtr, idx[0], idx[1], idx[2]);
      bFieldIn.setPtr(bFieldInPtr, idx[0], idx[1], idx[2]);
      temperatureIn.setPtr(temperatureInPtr, idx[0], idx[1], idx[2]);
      numDensityIn.setPtr(numDensityInPtr, idx[0], idx[1], idx[2]);

      // Compute scale factor to scale diffusion term so that total energy is zero
      double scaleFactor = -energyDragInPtr[0]/energyDiffInPtr[0];
      
      if (std::isinf(scaleFactor))
        continue;

      // Using scaleFactor, recompute CFL condition
      for (int configNode = 0; configNode < nlocal3d; configNode++)
      {
        // Keep track of max CFL number
        cfla = std::max( cfla, std::abs(4.0*alpha*numDensityInPtr[configNode]/(temperatureInPtr[configNode]*sqrt(temperatureInPtr[configNode]))
          *scaleFactor*temperatureInPtr[configNode]/speciesMass*dt/(grid.getDx(3)*grid.getDx(3))) );

        if (idx[4] < globalRgn.getUpper(4)-1)
        {
          double muCoord = cellCentroid[4] + 0.5*grid.getDx(4);
          double muTherm = scaleFactor*temperatureInPtr[configNode]/bFieldInPtr[configNode];
          cfla = std::max(cfla, 8.0*alpha*numDensityInPtr[configNode]/(temperatureInPtr[configNode]*sqrt(temperatureInPtr[configNode]))
            *muTherm*muCoord*dt/(grid.getDx(4)*grid.getDx(4)));
        }
      }

      distfDiff.setPtr(distfDiffPtr, idx);
      distfOut.setPtr(distfOutPtr, idx);
      // Set fOut to give the right energy at this node
      for (int nodeIndex = 0; nodeIndex < nlocal5d; nodeIndex++)
      {
        double startVal = distfOutPtr[nodeIndex];
        distfOutPtr[nodeIndex] = distfOutPtr[nodeIndex] + dt*scaleFactor*distfDiffPtr[nodeIndex];
      }
    }

    if (cfla > cflm)
      return Lucee::UpdaterStatus(false, dt*cfl/cfla);
    else
      return Lucee::UpdaterStatus(true, dt*cfl/cfla);
  }

  void
  SOLLenardBernsteinScaleCell5DUpdater::declareTypes()
  {
    // Input: Distribution function after diffusion step
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Input: Energy in each cell after drag step
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: Energy in each cell after drag step
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: Temperature in joules
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: Magnetic field
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: Dimensionally correct number density
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Distribution function with drag already applied and diffusion to be added on in update()
    this->appendOutVarType(typeid(Lucee::Field<5, double>));
  }

  void
  SOLLenardBernsteinScaleCell5DUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
  SOLLenardBernsteinScaleCell5DUpdater::evaluateFunction(Lucee::LuaState& L, double tm,
    std::vector<double>& res)
  {
    // push function object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, fnRef);
    // push variables on stack
    lua_pushnumber(L, tm);
    // call function
    if (lua_pcall(L, 1, res.size(), 0) != 0)
    {
      Lucee::Except lce("SOLLenardBernsteinScaleCell5DUpdater::evaluateFunction: ");
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
        throw Lucee::Except("SOLLenardBernsteinScaleCell5DUpdater::evaluateFunction: Return value not a number");
      res[res.size()+i] = lua_tonumber(L, i);
    }
    lua_pop(L, 1);
  }
}
