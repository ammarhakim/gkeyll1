/**
 * @file	LcSOLInitializeDensity.cpp
 *
 * @brief	Scale the distribution function at each node to have an exact density
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLInitializeDensity.h>
#include <LcGlobals.h>
// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  const char *SOLInitializeDensity::id = "SOLInitializeDensity";

  SOLInitializeDensity::SOLInitializeDensity()
    : fnRef(-1)
  {
  }

  void
  SOLInitializeDensity::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    // get function to evaluate
    fnRef = tbl.getFunctionRef("evaluate");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except("SOLInitializeDensity::readInput: Must specify element to use using 'basis5d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("SOLInitializeDensity::readInput: Must specify element to use using 'basis3d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis2d"))
      nodalBasis2d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis2d");
    else
      throw Lucee::Except("SOLInitializeDensity::readInput: Must specify element to use using 'basis2d'");

    if (tbl.hasNumber("scaleFactor"))
      scaleFactor = tbl.getNumber("scaleFactor");
    else scaleFactor = 1.0;
  }

  void
  SOLInitializeDensity::initialize()
  {
    UpdaterIfc::initialize();

    unsigned nlocal = nodalBasis5d->getNumNodes();
    std::vector<unsigned> zRef(nlocal), vRef(nlocal);

    unsigned nlocal3d = nodalBasis3d->getNumNodes();

    // Get a copy of the 5d nodal coordinates
    Lucee::Matrix<double> nodeCoordsLucee(nlocal, 5);
    nodalBasis5d->getNodalCoordinates(nodeCoordsLucee);
    Eigen::MatrixXd nodeCoords(nlocal, 5);
    copyLuceeToEigen(nodeCoordsLucee, nodeCoords);

    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // Used to figure out which nodes share the same location in configuration space
    double dxMin = grid.getDx(0);
    for (int d = 1; d < 3; d++)
      dxMin = std::min(dxMin, grid.getDx(d));

    // Find all nodes that share the same location as node zero in configuration space
    nodalStencil = std::vector<int>(nlocal);
    int stencilIndex = 0;
    for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
    {
      if (sameConfigCoords(0, nodeIndex, dxMin, nodeCoords) == true)
      {
        nodalStencil[stencilIndex] = nodeIndex;
        stencilIndex++;
      }
    }
    nodalStencil.resize(stencilIndex);

    // Compute matrix for v,mu integration
    int nlocal2d = nodalBasis2d->getNumNodes();
    Lucee::Matrix<double> tempMassMatrix(nlocal2d, nlocal2d);

    Eigen::MatrixXd integrationMatrix(nlocal2d, nlocal2d);
    nodalBasis2d->getMassMatrix(tempMassMatrix);
    copyLuceeToEigen(tempMassMatrix, integrationMatrix);
    integrationMatrix *= grid.getDx(3)*grid.getDx(4)/(grid.getDx(0)*grid.getDx(1));

    integrationVector = integrationMatrix.colwise().sum();
  }

  Lucee::UpdaterStatus
  SOLInitializeDensity::update(double t)
  {
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    const Lucee::Field<3, double>& bField = this->getInp<Lucee::Field<3, double> >(0);
    // Distribution function to be scaled
    Lucee::Field<5, double>& distf = this->getOut<Lucee::Field<5, double> >(0);
    
    Lucee::Region<5, int> globalRgn = grid.getGlobalRegion();
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();
    // get hold of Lua state object
    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;
    
    Lucee::ConstFieldPtr<double> bFieldPtr = bField.createConstPtr();
    Lucee::FieldPtr<double> distfPtr = distf.createPtr(); // Output pointer

    unsigned nlocal5d = nodalBasis5d->getNumNodes();
    unsigned nlocal3d = nodalBasis3d->getNumNodes();

    int idx[5];
    double xc[5];
    std::vector<double> result(1);
    Eigen::VectorXd distfReduced(nodalStencil.size());
    Lucee::Matrix<double> nodeCoords3dLucee(nlocal3d, 3);
    // Loop over each node in configuration space
    for (int ix = localRgn.getLower(0); ix < localRgn.getUpper(0); ix++)
    {
      for (int iy = localRgn.getLower(1); iy < localRgn.getUpper(1); iy++)
      {
        for (int iz = localRgn.getLower(2); iz < localRgn.getUpper(2); iz++)
        {
          idx[0] = ix;
          idx[1] = iy;
          idx[2] = iz;
          idx[3] = localRgn.getLower(3);
          idx[4] = localRgn.getLower(4);
          grid.setIndex(idx);
          grid.getCentroid(xc);
          bField.setPtr(bFieldPtr, idx[0], idx[1], idx[2]);

          // Get a copy of the 3d node coordinates
          nodalBasis3d->setIndex(idx[0], idx[1], idx[2]);
          nodalBasis3d->getNodalCoordinates(nodeCoords3dLucee);
          // Fill out absolute node list using local node coordinate list
          for (int configNode = 0; configNode < nlocal3d; configNode++)
          {
            // Get desired density at this configuration node
            evaluateFunction(*L, t, nodeCoords3dLucee, configNode, result);
            // Compute density of numerical solution at this configuration node
            double numericalDensity = 0.0;
            for (int iv = localRgn.getLower(3); iv < localRgn.getUpper(3); iv++)
            {
              idx[3] = iv;
              for (int imu = localRgn.getLower(4); imu < localRgn.getUpper(4); imu++)
              {
                idx[4] = imu;
                distf.setPtr(distfPtr, idx);
                for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                  distfReduced(nodeIndex) = distfPtr[nodalStencil[nodeIndex] + configNode];
                numericalDensity += integrationVector.dot(distfReduced);
              }
            }
            numericalDensity = numericalDensity*scaleFactor*bFieldPtr[configNode];
            //std::cout << "target = " << result[0] << ", num = " << numericalDensity << std::endl;

            // Loop through nodes a second time to scale to correct value
            for (int iv = localRgn.getLower(3); iv < localRgn.getUpper(3); iv++)
            {
              idx[3] = iv;
              for (int imu = localRgn.getLower(4); imu < localRgn.getUpper(4); imu++)
              {
                idx[4] = imu;
                distf.setPtr(distfPtr, idx);
                for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                  distfPtr[nodalStencil[nodeIndex] + configNode] = (result[0]/numericalDensity)*
                    distfPtr[nodalStencil[nodeIndex] + configNode];
              }
            }
          }
        }
      }
    }
   
    return Lucee::UpdaterStatus();
  }

  void
  SOLInitializeDensity::declareTypes()
  {
    // Magnetic field
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
    // Distribution function with correct density
    this->appendOutVarType(typeid(Lucee::Field<5, double>));
  }

  void
  SOLInitializeDensity::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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

  bool
  SOLInitializeDensity::sameConfigCoords(int srcIndex, int tarIndex, double dxMin,
    const Eigen::MatrixXd& nodeList)
  {
    for (int d = 0; d < 3; d++)
      if (std::fabs(nodeList(srcIndex,d)-nodeList(tarIndex,d)) > 1e-4*dxMin) 
        return false;
    return true;
  }

  void
  SOLInitializeDensity::evaluateFunction(Lucee::LuaState& L, double tm,
    const Lucee::Matrix<double> nodeList, unsigned nn, std::vector<double>& res)
  {
// push function object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, fnRef);
// push variables on stack
    for (unsigned i=0; i<NC; ++i)
      lua_pushnumber(L, nodeList(nn,i));
    lua_pushnumber(L, tm);
// call function
    if (lua_pcall(L, NC+1, res.size(), 0) != 0)
    {
      Lucee::Except lce("SOLInitializeDensity::evaluateFunction: ");
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
        throw Lucee::Except("SOLInitializeDensity::evaluateFunction: Return value not a number");
      res[res.size()+i] = lua_tonumber(L, i);
    }
    lua_pop(L, 1);
  }
}
