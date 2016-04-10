/**
 * @file	LcDistFuncReflectionBcUpdater.cpp
 *
 * @brief	Applies electrostatic logical sheath BCs to a 5D (electron) distribution function
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLogicalSheath5DUpdater.h>
#include <LcGlobals.h>
//#include <LcMathLib.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <limits>
#include <vector>

//etc includes
#include <quadrule.hpp>

namespace Lucee
{
  const char *LogicalSheath5DUpdater::id = "LogicalSheath5D";

  LogicalSheath5DUpdater::LogicalSheath5DUpdater()
  {
  }

  void
  LogicalSheath5DUpdater::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except("LogicalSheath5DUpdater::readInput: Must specify element to use using 'basis5d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("LogicalSheath5DUpdater::readInput: Must specify element to use using 'basis3d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis2d"))
      nodalBasis2d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis2d");
    else
      throw Lucee::Except("LogicalSheath5DUpdater::readInput: Must specify element to use using 'basis2d'");

    if (tbl.hasNumber("polyOrder"))
      polyOrder = tbl.getNumber("polyOrder");
    else
      throw Lucee::Except("LogicalSheath5DUpdater::readInput: Must specify basis function order using 'polyOrder'");
 
    // Factor to multiply all results by (like 2*pi*B/m to account v_perp -> mu integration
    if (tbl.hasNumber("scaleFactor"))
      scaleFactor = tbl.getNumber("scaleFactor");
    else scaleFactor = 1.0;

    // Mass of electrons
    if (tbl.hasNumber("elcMass"))
      elcMass = tbl.getNumber("elcMass");
    else
      throw Lucee::Except("LogicalSheath5DUpdater::readInput: Must specify electron mass using 'elcMass'");

    // Mass of ions
    if (tbl.hasNumber("ionMass"))
      ionMass = tbl.getNumber("ionMass");
    else
      throw Lucee::Except("LogicalSheath5DUpdater::readInput: Must specify ion mass using 'ionMass'");

    // Elementary charge
    if (tbl.hasNumber("eV"))
      eV = tbl.getNumber("eV");
    else
      throw Lucee::Except("LogicalSheath5DUpdater::readInput: Must specify electronvolt value using 'eV'");

    // Setting this to false will avoid a costly computation, but no sheath potential data will be available
    if (tbl.hasBool("computeCutoffVelocities"))
      computeCutoffVelocities = tbl.getBool("computeCutoffVelocities");
    else computeCutoffVelocities = true;
  }

  void
  LogicalSheath5DUpdater::initialize()
  {
    UpdaterIfc::initialize();

    unsigned nlocal = nodalBasis5d->getNumNodes();
    std::vector<unsigned> zRef(nlocal), vRef(nlocal);

    // Get reflection mapping after element has been reflected in z and v_para
    rotMap.resize(nlocal);
    nodalBasis5d->getUpperReflectingBcMapping(2, zRef);
    nodalBasis5d->getUpperReflectingBcMapping(3, vRef);
    for (int i = 0; i < nlocal; i++)
      rotMap[i] = vRef[zRef[i]];

    // Get a copy of the nodal coordinates
    Lucee::Matrix<double> nodeCoordsLucee(nlocal, 5);
    nodalBasis5d->getNodalCoordinates(nodeCoordsLucee);
    Eigen::MatrixXd nodeCoords(nlocal, 5);
    copyLuceeToEigen(nodeCoordsLucee, nodeCoords);

    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    double dxMin = grid.getDx(0);
    for (int d = 1; d < 3; d++)
      dxMin = std::min(dxMin, grid.getDx(d));

    // Find all nodes that share the same location as node zero
    // Will eventually need to do this at all nodes on lower surface
    // but can get away with doing this at one point for linear element
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

    // Construct a quadrature matrix to integrate over entire (v,mu) cell
    int integrationDegree = polyOrder + 1; // (standard linear basis function times one degree in v)
    unsigned numGaussPoints1d = (unsigned)((integrationDegree+1)/2.0 + 0.5);
    std::vector<double> gaussPoints1d(numGaussPoints1d);
    std::vector<double> gaussWeights1d(numGaussPoints1d);
    legendre_set(numGaussPoints1d, &gaussPoints1d[0], &gaussWeights1d[0]);

    int totalSurfQuadPoints = numGaussPoints1d*numGaussPoints1d;
    Eigen::MatrixXd gaussSurf(totalSurfQuadPoints, nodalStencil.size());
    gaussSurfCoords = Eigen::MatrixXd(totalSurfQuadPoints, 2);
    gaussSurfWeights = Eigen::VectorXd(totalSurfQuadPoints);
    double refCoord[5];
    refCoord[0] = -1;
    refCoord[1] = -1;
    refCoord[2] = -1;
    std::vector<double> basisAtPoint(nodalStencil.size());
    for (int gaussIndexOuter = 0; gaussIndexOuter < numGaussPoints1d; gaussIndexOuter++)
    {
      refCoord[3] = gaussPoints1d[gaussIndexOuter];
      for (int gaussIndexInner = 0; gaussIndexInner < numGaussPoints1d; gaussIndexInner++)
      {
        refCoord[4] = gaussPoints1d[gaussIndexInner];
        int linIndex = gaussIndexOuter*numGaussPoints1d + gaussIndexInner;
        // Store coordinate of quadrature point (on ref. element)
        gaussSurfCoords.row(linIndex) << refCoord[3], refCoord[4];
        // Store integration weight of quadrature point (real space)
        gaussSurfWeights(linIndex) = gaussWeights1d[gaussIndexOuter]*gaussWeights1d[gaussIndexInner];
        // Evaluate relevant basis functions at quadrature point
        nodalBasis5d->evalBasis(refCoord, basisAtPoint, nodalStencil);
        // Store results in matrix
        for (int nodeIndex = 0; nodeIndex < basisAtPoint.size(); nodeIndex++)
          gaussSurf(linIndex, nodeIndex) = basisAtPoint[nodeIndex];
      }
    }

    // Compute a matrix used to calculate the total parallel flux across a surface
    momentMatrix = Eigen::MatrixXd(nodalStencil.size(), nodalStencil.size());
    double weightScale = 0.5*grid.getDx(3)*0.5*grid.getDx(4);

    for (int i = 0; i < nodalStencil.size(); i++)
    {
      for (int j = 0; j < nodalStencil.size(); j++)
      {
        double integralResult = 0.0;
        for (int quadIndex = 0; quadIndex < totalSurfQuadPoints; quadIndex++)
          integralResult += weightScale*gaussSurfWeights(quadIndex)*
            gaussSurf(quadIndex, i)*gaussSurf(quadIndex, j);
        momentMatrix(i,j) = integralResult;
      }
    }

    // Store 5d mass matrix inverse
    Lucee::Matrix<double> tempMatrix(nlocal, nlocal);
    nodalBasis5d->getMassMatrix(tempMatrix);
    Eigen::MatrixXd massMatrix(nlocal, nlocal);
    copyLuceeToEigen(tempMatrix, massMatrix);
    // Store stiffness matrix
    nodalBasis5d->getGradStiffnessMatrix(3, tempMatrix);
    Eigen::MatrixXd gradStiffnessMatrix(nlocal, nlocal);
    copyLuceeToEigen(tempMatrix, gradStiffnessMatrix);
    // Compute and store differention matrix
    gradMatrix = massMatrix.inverse()*gradStiffnessMatrix.transpose();
    // Fill out the node numbers on lower and upper surfaces in z
    lowerEdgeNodeNums = std::vector<int>(nodalBasis3d->getNumSurfLowerNodes(2));
    upperEdgeNodeNums = std::vector<int>(nodalBasis3d->getNumSurfUpperNodes(2));
    nodalBasis3d->getSurfLowerNodeNums(2, lowerEdgeNodeNums);
    nodalBasis3d->getSurfUpperNodeNums(2, upperEdgeNodeNums);

    cutoffTolerance = 1e-4;
  }

  Lucee::UpdaterStatus
  LogicalSheath5DUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    const Lucee::Field<3, double>& ionFluxIn = this->getInp<Lucee::Field<3, double> >(0);
    const Lucee::Field<5, double>& hamilIn = this->getInp<Lucee::Field<5, double> >(1);
    // Output distribution function
    Lucee::Field<5, double>& distf = this->getOut<Lucee::Field<5, double> >(0);
    // Output sheath potential
    Lucee::Field<2, double>& phiSLower = this->getOut<Lucee::Field<2, double> >(1);
    Lucee::Field<2, double>& phiSUpper = this->getOut<Lucee::Field<2, double> >(2);

    Lucee::Region<5, int> globalRgn = grid.getGlobalRegion();
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();

    Lucee::ConstFieldPtr<double> ionFluxPtr = ionFluxIn.createConstPtr();
    Lucee::ConstFieldPtr<double> hamilPtr = hamilIn.createConstPtr();
    Lucee::FieldPtr<double> sknPtr = distf.createPtr(); // for skin-cell
    Lucee::FieldPtr<double> gstPtr = distf.createPtr(); // for ghost-cell
    Lucee::FieldPtr<double> phiSLowerPtr = phiSLower.createPtr(); // for lower sheath potential
    Lucee::FieldPtr<double> phiSUpperPtr = phiSUpper.createPtr(); // for upper sheath potential

    unsigned nlocal = nodalBasis5d->getNumNodes();
    unsigned nlocal3d = nodalBasis3d->getNumNodes();
    double cellCentroid[5];
    // index of skin cell
    int idx[5];
    // index of ghost cell
    int gstIdx[5];
    bool foundAllVc = true;
    const int maxIter = 50;

    // Need to find phiS on lower z plane if it is contained in the localRegion
    if (localRgn.getLower(2) == globalRgn.getLower(2))
    {
      // Outer loops are over the position cells
      for (int ix = localRgn.getLower(0); ix < localRgn.getUpper(0); ix++)
      {
        for (int iy = localRgn.getLower(1); iy < localRgn.getUpper(1); iy++)
        {
          idx[0] = ix;
          idx[1] = iy;
          idx[2] = localRgn.getLower(2);
          ionFluxIn.setPtr(ionFluxPtr, idx[0], idx[1], idx[2]);
          // Set output pointer for sheath potential
          phiSLower.setPtr(phiSLowerPtr, idx[0], idx[1]);
          
          gstIdx[0] = ix;
          gstIdx[1] = iy;
          gstIdx[2] = idx[2]-1;
          // Loop over all configuration space nodes contained in this element
          for (int configNode = 0; configNode < lowerEdgeNodeNums.size(); configNode++)
          {
            int configNodeIndex = lowerEdgeNodeNums[configNode];
            // This is the flux we want to match at this point
            double totalIonFluxAtNode = elcMass/ionMass*ionFluxPtr[configNodeIndex];
            double runningElcFluxAtNode = 0.0;
            bool foundCutoffCell = false;

            //if (totalIonFluxAtNode > 0.0)
            //  totalIonFluxAtNode = 0.0;

            // Need to loop over velocity space in this specific manner
            for (int ivSkin = localRgn.getLower(3), ivGhost = localRgn.getUpper(3)-1;
                ivSkin < localRgn.getUpper(3); ivSkin++, ivGhost--)
            {
              idx[3] = ivSkin;
              gstIdx[3] = ivGhost;

              if (foundCutoffCell == false)
              {
                // Need to check if this is the parallel velocity cutoff cell
                double elcFluxAtIv = 0.0;
                // At every velocity space index, need to compute flux over all cells in mu
                for (int iMu = localRgn.getLower(4); iMu < localRgn.getUpper(4); iMu++)
                {
                  idx[4] = iMu;
                  distf.setPtr(sknPtr, idx);
                  hamilIn.setPtr(hamilPtr, idx);
                  // First compute entire hamiltonian derivative
                  Eigen::VectorXd hamilFull(nlocal);
                  for (int i = 0; i < nlocal; i++)
                    hamilFull(i) = hamilPtr[i];
                  // Compute derivative of hamiltonian expressed in terms of basis functions
                  Eigen::VectorXd hamilDerivFull = gradMatrix*hamilFull;
                  // Get the coordinates of cell center
                  grid.setIndex(idx);
                  grid.getCentroid(cellCentroid);
                  // At this particular configuration space vertix, copy all
                  // nodes that occupy this location to a vector
                  Eigen::VectorXd distfReduced(nodalStencil.size());
                  Eigen::VectorXd hamilReduced(nodalStencil.size());
                  for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                  {
                    distfReduced(nodeIndex) = sknPtr[nodalStencil[nodeIndex] + configNodeIndex];
                    hamilReduced(nodeIndex) = hamilDerivFull(nodalStencil[nodeIndex] + configNodeIndex);
                  }
                  elcFluxAtIv += scaleFactor*distfReduced.dot(momentMatrix*hamilReduced);
                }
                // Check if we have exceed the total ion flux
                if (runningElcFluxAtNode + elcFluxAtIv < totalIonFluxAtNode)
                {
                  foundCutoffCell = true;
                  // Get the coordinates of cell center
                  idx[4] = localRgn.getLower(4);
                  grid.setIndex(idx);
                  grid.getCentroid(cellCentroid);
                  // Figure out fraction of cell that contains the excess flux
                  // (Flux over what is needed for equivalence with Gamma_i)
                  // exactResult is the the flux contribution integrated from lower edge of cell to v_c, a negative value
                  double exactResult = totalIonFluxAtNode - runningElcFluxAtNode; 
                  double excessFraction = (runningElcFluxAtNode + elcFluxAtIv - totalIonFluxAtNode)/elcFluxAtIv;
                  if (excessFraction < 0.0)
                    std::cout << "excessFraction negative" << std::endl;
                  // Search for the cutoff velocity ('b' here)
                  double a = -0.5*grid.getDx(3);
                  double b;
                  double relError;
                  int iterCount = 0;
                  std::vector<double> basisAtPoint(nodalStencil.size());
                  // Root bracket values
                  double lowerBound = -0.5*grid.getDx(3);
                  double upperBound = 0.5*grid.getDx(3);
                  // Gamma - Gamma_exact evaluated at bracketed values
                  double fl = 0.0-exactResult;
                  double fh = elcFluxAtIv-exactResult;

                  //std::cout << "elcFluxAtIv = " << elcFluxAtIv << std::endl;
                  //std::cout << "runningElcFluxAtNode = " << runningElcFluxAtNode << std::endl;
                  //std::cout << "totalIonFluxAtNode = " << totalIonFluxAtNode << std::endl;
                  // Ridders' Method (See Press 2007, page 453). Removed some of checks that
                  // the version in Numerical Recipes has because they appear to be unnecessary.
                  // Consider using Brent's method in the future (more complicated to implement).
                  /*
                  for (int iter = 0; iter < 60; iter++)
                  {
                    double xm = 0.5*(lowerBound + upperBound);
                    // Calculate function at xm
                    double fm = 0.0;
                    b = xm;
                    // Integration weight scale for this modified integral
                    double weightScale = 0.5*(b-a)*0.5*grid.getDx(4);
                    double refCoord[2];
                    
                    for (int gaussIndex = 0; gaussIndex < gaussSurfCoords.rows(); gaussIndex++)
                    {
                      refCoord[0] = ( 0.5*(b-a)*gaussSurfCoords(gaussIndex,0) + 0.5*(a+b) )/(0.5*grid.getDx(3));
                      refCoord[1] = gaussSurfCoords(gaussIndex,1);
                      nodalBasis2d->evalBasis(refCoord, basisAtPoint);

                      for (int iMu = localRgn.getLower(4); iMu < localRgn.getUpper(4); iMu++)
                      {
                        idx[4] = iMu;
                        distf.setPtr(sknPtr, idx);
                        // Get the coordinates of cell center
                        grid.setIndex(idx);
                        grid.getCentroid(cellCentroid);
                        // Compute distribution function at quadrature point
                        double fAtPoint = 0.0;
                        for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                          fAtPoint += sknPtr[nodalStencil[nodeIndex] + configNodeIndex]*basisAtPoint[nodeIndex];

                        fm += weightScale*gaussSurfWeights(gaussIndex)*
                          (cellCentroid[3] + refCoord[0]*0.5*grid.getDx(3))*fAtPoint;
                      }
                    }

                    fm -= exactResult;

                    double s = sqrt(fm*fm - fl*fh);
                    // Updating formula
                    double xNew = xm + (xm-lowerBound)*( (fl >= fh ? 1.0 : -1.0)*fm/s );
                    // Evaluate f at xNew
                    double fNew = 0.0;

                    b = xNew;
                    // Integration weight scale for this modified integral
                    weightScale = 0.5*(b-a)*0.5*grid.getDx(4);
                    
                    for (int gaussIndex = 0; gaussIndex < gaussSurfCoords.rows(); gaussIndex++)
                    {
                      refCoord[0] = ( 0.5*(b-a)*gaussSurfCoords(gaussIndex,0) + 0.5*(a+b) )/(0.5*grid.getDx(3));
                      refCoord[1] = gaussSurfCoords(gaussIndex,1);
                      nodalBasis2d->evalBasis(refCoord, basisAtPoint);

                      for (int iMu = localRgn.getLower(4); iMu < localRgn.getUpper(4); iMu++)
                      {
                        idx[4] = iMu;
                        distf.setPtr(sknPtr, idx);
                        // Get the coordinates of cell center
                        grid.setIndex(idx);
                        grid.getCentroid(cellCentroid);
                        // Compute distribution function at quadrature point
                        double fAtPoint = 0.0;
                        for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                          fAtPoint += sknPtr[nodalStencil[nodeIndex] + configNodeIndex]*basisAtPoint[nodeIndex];

                        fNew += weightScale*gaussSurfWeights(gaussIndex)*
                          (cellCentroid[3] + refCoord[0]*0.5*grid.getDx(3))*fAtPoint;
                      }
                    }

                    fNew-=exactResult;

                    // fNew already is the difference between Gamma_tar and Gamma_exact
                    relError = fNew/exactResult;

                    
                    if (fabs(relError) < cutoffTolerance )
                    {
                      //std::cout << "lFinal relError = " << relError << std::endl;
                      //std::cout << "lFinal b = " << cellCentroid[3] + b << std::endl;
                      //std::cout << "lexactResult = " << exactResult << std::endl;
                      //std::cout << "literCount = " << iter << std::endl;
                      //std::cout << "lExcessFraction = " << excessFraction << std::endl;
                      break;
                    }

                    // Update root brackets, making use of monotonicity of function
                    if (relError > 0)
                    {
                      upperBound = b;
                      fh = fNew;
                    }
                    else
                    {
                      lowerBound = b;
                      fl = fNew;
                    }

                    if (iter > 20)
                    {
                      //std::cout << "iter = " << iter << std::endl;
                      //std::cout << "fl = " << fl << std::endl;
                      //std::cout << "fh = " << fh << std::endl;
                      //std::cout << "relError = " << relError << std::endl;
                      foundAllVc = false;
                    }
                  }*/

                  //if (exactResult == 0.0)
                  //  b = -0.5*grid.getDx(3);
                  //else
                  {
                    do
                    {
                      b = 0.5*(upperBound + lowerBound);
                      // Integration weight scale for this modified integral
                      double weightScale = 0.5*(b-a)*0.5*grid.getDx(4);
                      double integralResult = 0.0;
                      double refCoord[2];
                      
                      for (int gaussIndex = 0; gaussIndex < gaussSurfCoords.rows(); gaussIndex++)
                      {
                        refCoord[0] = ( 0.5*(b-a)*gaussSurfCoords(gaussIndex,0) + 0.5*(a+b) )/(0.5*grid.getDx(3));
                        refCoord[1] = gaussSurfCoords(gaussIndex,1);
                        nodalBasis2d->evalBasis(refCoord, basisAtPoint);

                        for (int iMu = localRgn.getLower(4); iMu < localRgn.getUpper(4); iMu++)
                        {
                          idx[4] = iMu;
                          distf.setPtr(sknPtr, idx);
                          hamilIn.setPtr(hamilPtr, idx);
                          // First compute entire hamiltonian derivative
                          Eigen::VectorXd hamilFull(nlocal);
                          for (int i = 0; i < nlocal; i++)
                            hamilFull(i) = hamilPtr[i];
                          // Compute derivative of hamiltonian expressed in terms of basis functions
                          Eigen::VectorXd hamilDerivFull = gradMatrix*hamilFull;
                          // Get the coordinates of cell center
                          grid.setIndex(idx);
                          grid.getCentroid(cellCentroid);
                          // Compute distribution function at quadrature point
                          double fAtPoint = 0.0;
                          double gradHAtPoint = 0.0;
                          for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                          {
                            fAtPoint += sknPtr[nodalStencil[nodeIndex] + configNodeIndex]*basisAtPoint[nodeIndex];
                            gradHAtPoint += hamilDerivFull(nodalStencil[nodeIndex] + configNodeIndex)*basisAtPoint[nodeIndex];
                          }

                          integralResult += weightScale*gaussSurfWeights(gaussIndex)*
                            scaleFactor*gradHAtPoint*fAtPoint;
                        }
                      }

                      // integralResult will be more negative than exactResult if guess
                      // is more "left" than true answer
                      relError = (integralResult - exactResult)/exactResult;

                      if (relError > 0)
                        upperBound = b;
                      else
                        lowerBound = b;

                      iterCount++;

                    } while ( fabs(relError) > cutoffTolerance && iterCount < maxIter);

                    if (iterCount == maxIter)
                    {
                      foundAllVc = false;
                      // do the search again for debug purposes, printing out more information
                      lowerBound = -0.5*grid.getDx(3);
                      upperBound = 0.5*grid.getDx(3);
                      iterCount = 0;
                      do
                      {
                        b = 0.5*(upperBound + lowerBound);
                        // Integration weight scale for this modified integral
                        double weightScale = 0.5*(b-a)*0.5*grid.getDx(4);
                        double integralResult = 0.0;
                        double refCoord[2];
                        
                        for (int gaussIndex = 0; gaussIndex < gaussSurfCoords.rows(); gaussIndex++)
                        {
                          refCoord[0] = ( 0.5*(b-a)*gaussSurfCoords(gaussIndex,0) + 0.5*(a+b) )/(0.5*grid.getDx(3));
                          refCoord[1] = gaussSurfCoords(gaussIndex,1);
                          nodalBasis2d->evalBasis(refCoord, basisAtPoint);

                          for (int iMu = localRgn.getLower(4); iMu < localRgn.getUpper(4); iMu++)
                          {
                            idx[4] = iMu;
                            distf.setPtr(sknPtr, idx);
                            hamilIn.setPtr(hamilPtr, idx);
                            // First compute entire hamiltonian derivative
                            Eigen::VectorXd hamilFull(nlocal);
                            for (int i = 0; i < nlocal; i++)
                              hamilFull(i) = hamilPtr[i];
                            // Compute derivative of hamiltonian expressed in terms of basis functions
                            Eigen::VectorXd hamilDerivFull = gradMatrix*hamilFull;
                            // Get the coordinates of cell center
                            grid.setIndex(idx);
                            grid.getCentroid(cellCentroid);
                            // Compute distribution function at quadrature point
                            double fAtPoint = 0.0;
                            double gradHAtPoint = 0.0;
                            for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                            {
                              fAtPoint += sknPtr[nodalStencil[nodeIndex] + configNodeIndex]*basisAtPoint[nodeIndex];
                              gradHAtPoint += hamilDerivFull(nodalStencil[nodeIndex] + configNodeIndex)*basisAtPoint[nodeIndex];
                            }
                            std::cout << "fAtPoint = " << fAtPoint << std::endl;

                            integralResult += weightScale*gaussSurfWeights(gaussIndex)*
                              scaleFactor*gradHAtPoint*fAtPoint;
                          }
                        }
                        // integralResult will be more negative than exactResult if guess
                        // is more "left" than true answer
                        relError = (integralResult - exactResult)/exactResult;

                        std::cout << "lower iter = " << iterCount << std::endl;
                        std::cout << "idx = " << idx[0] << "," << idx[1] << "," << idx[2] << "," << idx[3] << std::endl;
                        std::cout << "relError = " << relError << std::endl;
                        std::cout << "xm = " << b << std::endl;
                        std::cout << "integralResult = " << integralResult << std::endl;
                        std::cout << "exactResult = " << exactResult << std::endl;
                        std::cout << "lowerBound = " << lowerBound << std::endl;
                        std::cout << "upperBound = " << upperBound << std::endl;
                        std::cout << "totalIonFluxAtNode = " << totalIonFluxAtNode << std::endl;
                        std::cout << "runningElcFluxAtNode = " << runningElcFluxAtNode << std::endl << std::endl;

                        if (relError > 0)
                          upperBound = b;
                        else
                          lowerBound = b;

                        iterCount++;

                      } while ( fabs(relError) > cutoffTolerance && iterCount < maxIter);
                    }
                  }

                  // Store result in the appropriate 2d field
                  phiSLowerPtr[configNode] = 0.5*elcMass*(cellCentroid[3] + b)*(cellCentroid[3] + b)/eV;
                  //std::cout << "idx = " << idx[0] << "," << idx[1] << "," << idx[2] << "," << idx[3] << std::endl;
                  //std::cout << "phiSLowerPtr[" << configNode << "] = " << phiSLowerPtr[configNode] << std::endl;

                  // Scale (only for cutoff cells) and reflect distribution function by copying into ghost cell
                  for (int iMu = localRgn.getLower(4); iMu < localRgn.getUpper(4); iMu++)
                  {
                    idx[4] = iMu;
                    distf.setPtr(sknPtr, idx);
                    gstIdx[4] = iMu;
                    distf.setPtr(gstPtr, gstIdx);
                    // It's important not to reflect every node in the cell!
                    for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                    {
                      // rotMap[k] tells you what sknPtr node to get data to put into 'k'
                      // Equivalently, it tells you what gstPtr node to put data from sknPtr node 'k'
                      int sknNode = nodalStencil[nodeIndex] + configNodeIndex;
                      int gstNode = rotMap[sknNode];
                      gstPtr[gstNode] = excessFraction*sknPtr[sknNode];
                    }
                  }
                }
                else 
                {
                  runningElcFluxAtNode += elcFluxAtIv;
                  // Make sure electron distribution function in ghost cell is zeroed out
                  for (int iMu = localRgn.getLower(4); iMu < localRgn.getUpper(4); iMu++)
                  {
                    idx[4] = iMu;
                    distf.setPtr(sknPtr, idx);
                    gstIdx[4] = iMu;
                    distf.setPtr(gstPtr, gstIdx);
                    for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                    {
                      // rotMap[k] tells you what sknPtr node to get data to put into 'k'
                      // Equivalently, it tells you what gstPtr node to put data from sknPtr node 'k'
                      int sknNode = nodalStencil[nodeIndex] + configNodeIndex;
                      int gstNode = rotMap[sknNode];
                      gstPtr[gstNode] = 0.0;
                    }
                  }
                }
              }
              else
              {
                // Cutoff cell has already been found
                idx[4] = localRgn.getLower(4);
                // Get the coordinates of cell center
                grid.setIndex(idx);
                grid.getCentroid(cellCentroid);
                // Check to see if we need to do any copying and reflecting
                if (cellCentroid[3] < 0.0)
                {
                  // We've already found the cutoff velocity at this config node.
                  // Reflect distribution function by copying into ghost cell
                  for (int iMu = localRgn.getLower(4); iMu < localRgn.getUpper(4); iMu++)
                  {
                    idx[4] = iMu;
                    distf.setPtr(sknPtr, idx);
                    gstIdx[4] = iMu;
                    distf.setPtr(gstPtr, gstIdx);
                    // It's important not to reflect every node in the cell! Just the ones
                    // at this particular configNode
                    for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                    {
                      // rotMap[k] tells you what sknPtr node to get data to put into 'k'
                      // Equivalently, it tells you what gstPtr node to put data from sknPtr node 'k'
                      int sknNode = nodalStencil[nodeIndex] + configNodeIndex;
                      int gstNode = rotMap[sknNode];
                      gstPtr[gstNode] = sknPtr[sknNode];
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    // The following code is essentially the same as the above, with some subtle changes
    // in indexing and conditionals.
    // Need to find phiS on upper z plane if it is contained in the localRegion
    if (localRgn.getUpper(2) == globalRgn.getUpper(2))
    {
      // Outer loops are over the position cells
      for (int ix = localRgn.getLower(0); ix < localRgn.getUpper(0); ix++)
      {
        for (int iy = localRgn.getLower(1); iy < localRgn.getUpper(1); iy++)
        {
          idx[0] = ix;
          idx[1] = iy;
          idx[2] = localRgn.getUpper(2)-1;
          ionFluxIn.setPtr(ionFluxPtr, idx[0], idx[1], idx[2]);
          // Set output pointer for sheath potential
          phiSUpper.setPtr(phiSUpperPtr, idx[0], idx[1]);
          
          gstIdx[0] = ix;
          gstIdx[1] = iy;
          gstIdx[2] = idx[2]+1;
          // Loop over all configuration space nodes contained in this element
          for (int configNode = 0; configNode < upperEdgeNodeNums.size(); configNode++)
          {
            int configNodeIndex = upperEdgeNodeNums[configNode];
            // This is the flux we want to match at this point
            double totalIonFluxAtNode = elcMass/ionMass*ionFluxPtr[configNodeIndex];
            double runningElcFluxAtNode = 0.0;
            bool foundCutoffCell = false;

            // Need to loop over velocity space in this specific manner
            for (int ivSkin = localRgn.getUpper(3)-1, ivGhost = localRgn.getLower(3);
                ivSkin >= localRgn.getLower(3); ivSkin--, ivGhost++)
            {
              idx[3] = ivSkin;
              gstIdx[3] = ivGhost;

              if (foundCutoffCell == false)
              {
                // Need to check if this is the parallel velocity cutoff cell
                double elcFluxAtIv = 0.0;
                // At every velocity space index, need to compute flux over all cells in mu
                for (int iMu = localRgn.getLower(4); iMu < localRgn.getUpper(4); iMu++)
                {
                  idx[4] = iMu;
                  distf.setPtr(sknPtr, idx);
                  hamilIn.setPtr(hamilPtr, idx);
                  // First compute entire hamiltonian derivative
                  Eigen::VectorXd hamilFull(nlocal);
                  for (int i = 0; i < nlocal; i++)
                    hamilFull(i) = hamilPtr[i];
                  // Compute derivative of hamiltonian expressed in terms of basis functions
                  Eigen::VectorXd hamilDerivFull = gradMatrix*hamilFull;
                  Eigen::VectorXd hamilReduced(nodalStencil.size());
                  // Get the coordinates of cell center
                  grid.setIndex(idx);
                  grid.getCentroid(cellCentroid);
                  // At this particular configuration space vertix, copy all
                  // nodes that occupy this location to a vector
                  Eigen::VectorXd distfReduced(nodalStencil.size());
                  for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                  {
                    distfReduced(nodeIndex) = sknPtr[nodalStencil[nodeIndex] + configNodeIndex];
                    hamilReduced(nodeIndex) = hamilDerivFull(nodalStencil[nodeIndex] + configNodeIndex);
                  }
                  elcFluxAtIv += scaleFactor*distfReduced.dot(momentMatrix*hamilReduced);
                }
                // Check if we have exceed the total ion flux
                if (runningElcFluxAtNode + elcFluxAtIv > totalIonFluxAtNode)
                {
                  foundCutoffCell = true;
                  // Get the coordinates of cell center
                  idx[4] = localRgn.getLower(4);
                  grid.setIndex(idx);
                  grid.getCentroid(cellCentroid);
                  // (Flux over what is needed for equivalence with Gamma_i)
                  double exactResult = totalIonFluxAtNode - runningElcFluxAtNode; // the flux contribution up to v_c
                  // Figure out fraction of cell is excess, above what is needed for equivalent
                  double excessFraction = (runningElcFluxAtNode + elcFluxAtIv - totalIonFluxAtNode)/elcFluxAtIv;
                  if (excessFraction < 0.0)
                    std::cout << "excessFraction negative" << std::endl;
                  // Search for the cutoff velocity ('a' here)
                  double a = -0.5*grid.getDx(3);
                  double b = 0.5*grid.getDx(3);
                  double relError;
                  int iterCount = 0;
                  std::vector<double> basisAtPoint(nodalStencil.size());
                  // Root bracket values
                  double lowerBound = -0.5*grid.getDx(3);
                  double upperBound = 0.5*grid.getDx(3);
                  // Gamma - Gamma_exact evaluated at bracketed values
                  // Difference between flux contribution from v_c vs exact result for this cell
                  double fl = elcFluxAtIv-exactResult;
                  double fh = -exactResult;

                  //std::cout << "idx = " << idx[0] << "," << idx[1] << "," << idx[2] << "," << idx[3] << std::endl;
                  //std::cout << "fl = " << fl << std::endl; 
                  //std::cout << "fh = " << fh << std::endl; 
                  //std::cout << "totalIonFluxAtNode = " << totalIonFluxAtNode << std::endl;

                  // Ridders' Method (See Press 2007, page 453). Removed some of checks that
                  // the version in Numerical Recipes has because they appear to be unnecessary.
                  // Consider using Brent's method in the future (more complicated to implement).
                  /*
                  for (int iter = 0; iter < 60; iter++)
                  {
                    double xm = 0.5*(lowerBound + upperBound);
                    // Calculate function at xm
                    double fm = 0.0;
                    a = xm;
                    // Integration weight scale for this modified integral
                    double weightScale = 0.5*(b-a)*0.5*grid.getDx(4);
                    double refCoord[2];

                    for (int gaussIndex = 0; gaussIndex < gaussSurfCoords.rows(); gaussIndex++)
                    {
                      refCoord[0] = ( 0.5*(b-a)*gaussSurfCoords(gaussIndex,0) + 0.5*(a+b) )/(0.5*grid.getDx(3));
                      refCoord[1] = gaussSurfCoords(gaussIndex,1);
                      nodalBasis2d->evalBasis(refCoord, basisAtPoint);

                      for (int iMu = localRgn.getLower(4); iMu < localRgn.getUpper(4); iMu++)
                      {
                        idx[4] = iMu;
                        distf.setPtr(sknPtr, idx);
                        // Get the coordinates of cell center
                        grid.setIndex(idx);
                        grid.getCentroid(cellCentroid);
                        // Compute distribution function at quadrature point
                        double fAtPoint = 0.0;
                        for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                          fAtPoint += sknPtr[nodalStencil[nodeIndex] + configNodeIndex]*basisAtPoint[nodeIndex];

                        fm += weightScale*gaussSurfWeights(gaussIndex)*
                          (cellCentroid[3] + refCoord[0]*0.5*grid.getDx(3))*fAtPoint;
                      }
                    }

                    fm -= exactResult;

                    double s = sqrt(fm*fm - fl*fh);
                    // Updating formula
                    double xNew = xm + (xm-lowerBound)*( (fl >= fh ? 1.0 : -1.0)*fm/s );
                    // Evaluate f at xNew
                    double fNew = 0.0;

                    a = xNew;
                    // Integration weight scale for this modified integral
                    weightScale = 0.5*(b-a)*0.5*grid.getDx(4);
                    
                    for (int gaussIndex = 0; gaussIndex < gaussSurfCoords.rows(); gaussIndex++)
                    {
                      refCoord[0] = ( 0.5*(b-a)*gaussSurfCoords(gaussIndex,0) + 0.5*(a+b) )/(0.5*grid.getDx(3));
                      refCoord[1] = gaussSurfCoords(gaussIndex,1);
                      nodalBasis2d->evalBasis(refCoord, basisAtPoint);

                      for (int iMu = localRgn.getLower(4); iMu < localRgn.getUpper(4); iMu++)
                      {
                        idx[4] = iMu;
                        distf.setPtr(sknPtr, idx);
                        // Get the coordinates of cell center
                        grid.setIndex(idx);
                        grid.getCentroid(cellCentroid);
                        // Compute distribution function at quadrature point
                        double fAtPoint = 0.0;
                        for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                          fAtPoint += sknPtr[nodalStencil[nodeIndex] + configNodeIndex]*basisAtPoint[nodeIndex];

                        fNew += weightScale*gaussSurfWeights(gaussIndex)*
                          (cellCentroid[3] + refCoord[0]*0.5*grid.getDx(3))*fAtPoint;
                      }
                    }

                    fNew-=exactResult;

                    // fNew already is the difference between Gamma_tar and Gamma_exact
                    relError = fNew/exactResult;

                    
                    if (fabs(relError) < cutoffTolerance )
                    {
                      //std::cout << "uFinal relError = " << relError << std::endl;
                      //std::cout << "uFinal b = " << cellCentroid[3] + a << std::endl;
                      //std::cout << "uexactResult = " << exactResult << std::endl;
                      //std::cout << "uiterCount = " << iter << std::endl;
                      break;
                    }

                    // Update root brackets, making use of monotonicity of function
                    if (relError < 0)
                    {
                      upperBound = a;
                      fh = fNew;
                    }
                    else
                    {
                      lowerBound = a;
                      fl = fNew;
                    }

                    std::cout << "upper iter = " << iter << std::endl;
                    std::cout << "upper fl = " << fl << std::endl;
                    std::cout << "upper lowerBound = " << lowerBound << std::endl;
                    std::cout << "upper fh = " << fh << std::endl;
                    std::cout << "upper upperBound = " << upperBound << std::endl;
                    std::cout << "upper relError = " << relError << std::endl;
                    // check to see if we are at the last iteration
                    if (iter > 20)
                    {
                      foundAllVc = false;
                    }
                  }*/

                  // This is a bisection search for the exact 'a', keeping 'b' fixed
                  //if (exactResult == 0.0)
                  //{
                  //  a = 0.5*grid.getDx(3);
                  //}
                  //else
                  {
                    do
                    {
                      a = 0.5*(upperBound + lowerBound);
                      // Integration weight scale for this modified integral
                      double weightScale = 0.5*(b-a)*0.5*grid.getDx(4);
                      double integralResult = 0.0;
                      double refCoord[2];
                      
                      for (int gaussIndex = 0; gaussIndex < gaussSurfCoords.rows(); gaussIndex++)
                      {
                        refCoord[0] = ( 0.5*(b-a)*gaussSurfCoords(gaussIndex,0) + 0.5*(a+b) )/(0.5*grid.getDx(3));
                        refCoord[1] = gaussSurfCoords(gaussIndex,1);
                        nodalBasis2d->evalBasis(refCoord, basisAtPoint);

                        for (int iMu = localRgn.getLower(4); iMu < localRgn.getUpper(4); iMu++)
                        {
                          idx[4] = iMu;
                          distf.setPtr(sknPtr, idx);
                          hamilIn.setPtr(hamilPtr, idx);
                          // First compute entire hamiltonian derivative
                          Eigen::VectorXd hamilFull(nlocal);
                          for (int i = 0; i < nlocal; i++)
                            hamilFull(i) = hamilPtr[i];
                          // Compute derivative of hamiltonian expressed in terms of basis functions
                          Eigen::VectorXd hamilDerivFull = gradMatrix*hamilFull;
                          // Get the coordinates of cell center
                          grid.setIndex(idx);
                          grid.getCentroid(cellCentroid);
                          // Compute distribution function at quadrature point
                          double fAtPoint = 0.0;
                          double gradHAtPoint = 0.0;
                          for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                          {
                            fAtPoint += sknPtr[nodalStencil[nodeIndex] + configNodeIndex]*basisAtPoint[nodeIndex];
                            gradHAtPoint += hamilDerivFull(nodalStencil[nodeIndex] + configNodeIndex)*basisAtPoint[nodeIndex];
                          }

                          integralResult += weightScale*gaussSurfWeights(gaussIndex)*
                            scaleFactor*gradHAtPoint*fAtPoint;
                        }
                      }

                      relError = (integralResult - exactResult)/exactResult;

                      if (relError > 0)
                        lowerBound = a;
                      else
                        upperBound = a;

                      iterCount++;

                    } while ( fabs(relError) > cutoffTolerance && iterCount < maxIter);

                    if (iterCount == maxIter)
                    {
                      foundAllVc = false;
                      // do the search again for debug purposes, printing out more information
                      lowerBound = -0.5*grid.getDx(3);
                      upperBound = 0.5*grid.getDx(3);
                      iterCount = 0;

                      do
                      {
                        a = 0.5*(upperBound + lowerBound);
                        // Integration weight scale for this modified integral
                        double weightScale = 0.5*(b-a)*0.5*grid.getDx(4);
                        double integralResult = 0.0;
                        double refCoord[2];
                        
                        for (int gaussIndex = 0; gaussIndex < gaussSurfCoords.rows(); gaussIndex++)
                        {
                          refCoord[0] = ( 0.5*(b-a)*gaussSurfCoords(gaussIndex,0) + 0.5*(a+b) )/(0.5*grid.getDx(3));
                          refCoord[1] = gaussSurfCoords(gaussIndex,1);
                          nodalBasis2d->evalBasis(refCoord, basisAtPoint);

                          for (int iMu = localRgn.getLower(4); iMu < localRgn.getUpper(4); iMu++)
                          {
                            idx[4] = iMu;
                            distf.setPtr(sknPtr, idx);
                            hamilIn.setPtr(hamilPtr, idx);
                            // First compute entire hamiltonian derivative
                            Eigen::VectorXd hamilFull(nlocal);
                            for (int i = 0; i < nlocal; i++)
                              hamilFull(i) = hamilPtr[i];
                            // Compute derivative of hamiltonian expressed in terms of basis functions
                            Eigen::VectorXd hamilDerivFull = gradMatrix*hamilFull;
                            // Get the coordinates of cell center
                            grid.setIndex(idx);
                            grid.getCentroid(cellCentroid);
                            // Compute distribution function at quadrature point
                            double fAtPoint = 0.0;
                            double gradHAtPoint = 0.0;
                            for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                            {
                              fAtPoint += sknPtr[nodalStencil[nodeIndex] + configNodeIndex]*basisAtPoint[nodeIndex];
                              gradHAtPoint += hamilDerivFull(nodalStencil[nodeIndex] + configNodeIndex)*basisAtPoint[nodeIndex];
                            }
                            std::cout << "fAtPoint = " << fAtPoint << std::endl;
                            integralResult += weightScale*gaussSurfWeights(gaussIndex)*
                              scaleFactor*gradHAtPoint*fAtPoint;
                          }
                        }

                        relError = (integralResult - exactResult)/exactResult;

                        std::cout << "upper iter = " << iterCount << std::endl;
                        std::cout << "idx = " << idx[0] << "," << idx[1] << "," << idx[2] << "," << idx[3] << std::endl;
                        std::cout << "relError = " << relError << std::endl;
                        std::cout << "xm = " << a << std::endl;
                        std::cout << "integralResult = " << integralResult << std::endl;
                        std::cout << "exactResult = " << exactResult << std::endl;
                        std::cout << "lowerBound = " << lowerBound << std::endl;
                        std::cout << "upperBound = " << upperBound << std::endl;
                        std::cout << "totalIonFluxAtNode = " << totalIonFluxAtNode << std::endl;
                        std::cout << "runningElcFluxAtNode = " << runningElcFluxAtNode << std::endl << std::endl;

                        if (relError > 0)
                          lowerBound = a;
                        else
                          upperBound = a;

                        iterCount++;

                      } while ( fabs(relError) > cutoffTolerance && iterCount < maxIter);
                    }
                  }

                  // Store result in the appropriate 2d field
                  phiSUpperPtr[configNode] = 0.5*elcMass*(cellCentroid[3] + a)*(cellCentroid[3] + a)/eV;

                  // Scale (only for cutoff cells) and reflect distribution function by copying into ghost cell
                  for (int iMu = localRgn.getLower(4); iMu < localRgn.getUpper(4); iMu++)
                  {
                    idx[4] = iMu;
                    distf.setPtr(sknPtr, idx);
                    gstIdx[4] = iMu;
                    distf.setPtr(gstPtr, gstIdx);
                    // It's important not to reflect every node in the cell!
                    for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                    {
                      // rotMap[k] tells you what sknPtr node to get data to put into 'k'
                      // Equivalently, it tells you what gstPtr node to put data from sknPtr node 'k'
                      int sknNode = nodalStencil[nodeIndex] + configNodeIndex;
                      int gstNode = rotMap[sknNode];
                      gstPtr[gstNode] = excessFraction*sknPtr[sknNode];
                    }
                  }
                }
                else
                {
                  runningElcFluxAtNode += elcFluxAtIv;
                  // Make sure electron distribution function in ghost cell is zeroed out
                  for (int iMu = localRgn.getLower(4); iMu < localRgn.getUpper(4); iMu++)
                  {
                    idx[4] = iMu;
                    distf.setPtr(sknPtr, idx);
                    gstIdx[4] = iMu;
                    distf.setPtr(gstPtr, gstIdx);
                    for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                    {
                      // rotMap[k] tells you what sknPtr node to get data to put into 'k'
                      // Equivalently, it tells you what gstPtr node to put data from sknPtr node 'k'
                      int sknNode = nodalStencil[nodeIndex] + configNodeIndex;
                      int gstNode = rotMap[sknNode];
                      gstPtr[gstNode] = 0.0;
                    }
                  }
                }
              }
              else
              {
                idx[4] = localRgn.getLower(4);
                // Get the coordinates of cell center
                grid.setIndex(idx);
                grid.getCentroid(cellCentroid);
                // Check to see if we need to do any copying and reflecting
                if (cellCentroid[3] > 0.0)
                {
                  // We've already found the cutoff velocity at this config node.
                  // Reflect distribution function by copying into ghost cell
                  for (int iMu = localRgn.getLower(4); iMu < localRgn.getUpper(4); iMu++)
                  {
                    idx[4] = iMu;
                    distf.setPtr(sknPtr, idx);
                    gstIdx[4] = iMu;
                    distf.setPtr(gstPtr, gstIdx);
                    // It's important not to reflect every node in the cell! Just the ones
                    // at this particular configNode
                    for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                    {
                      // rotMap[k] tells you what sknPtr node to get data to put into 'k'
                      // Equivalently, it tells you what gstPtr node to put data from sknPtr node 'k'
                      int sknNode = nodalStencil[nodeIndex] + configNodeIndex;
                      int gstNode = rotMap[sknNode];
                      gstPtr[gstNode] = sknPtr[sknNode];
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    //return foundAllVc;
    if (foundAllVc == false)
      return Lucee::UpdaterStatus(false, 0.0);
    else return Lucee::UpdaterStatus();
  }

  void
  LogicalSheath5DUpdater::declareTypes()
  {
    // Moments of the ion distribution function on both edges
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Electron hamiltonian
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Distribution function output (electrons)
    this->appendOutVarType(typeid(Lucee::Field<5, double>));
    // Output sheath potential (lower z surface)
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
    // Output sheath potential (upper z surface)
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }

  void
  LogicalSheath5DUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
  LogicalSheath5DUpdater::sameConfigCoords(int srcIndex, int tarIndex, double dxMin,
    const Eigen::MatrixXd& nodeList)
  {
    for (int d = 0; d < 3; d++)
      if (std::fabs(nodeList(srcIndex,d)-nodeList(tarIndex,d)) > 1e-4*dxMin) 
        return false;
    return true;
  }
}
