/**
 * @file	LcDistFuncReflectionBcUpdater.cpp
 *
 * @brief	Applies electrostatic logical sheath BCs to a 5D distribution function
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBiasedSheath5DUpdater.h>

//etc includes
#include <quadrule.hpp>

namespace Lucee
{
  const char *BiasedSheath5DUpdater::id = "BiasedSheath5D";

  BiasedSheath5DUpdater::BiasedSheath5DUpdater()
  {
  }

  void
  BiasedSheath5DUpdater::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except("BiasedSheath5DUpdater::readInput: Must specify element to use using 'basis5d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("BiasedSheath5DUpdater::readInput: Must specify element to use using 'basis3d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis2d"))
      nodalBasis2d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis2d");
    else
      throw Lucee::Except("BiasedSheath5DUpdater::readInput: Must specify element to use using 'basis2d'");

    if (tbl.hasNumber("polyOrder"))
      polyOrder = tbl.getNumber("polyOrder");
    else
      throw Lucee::Except("BiasedSheath5DUpdater::readInput: Must specify basis function order using 'polyOrder'");
 
    // Factor to multiply all results by (like 2*pi*B/m to account v_perp -> mu integration
    if (tbl.hasNumber("scaleFactor"))
      scaleFactor = tbl.getNumber("scaleFactor");
    else scaleFactor = 1.0;

    // Mass of species
    if (tbl.hasNumber("speciesMass"))
      speciesMass = tbl.getNumber("speciesMass");
    else
      throw Lucee::Except("BiasedSheath5DUpdater::readInput: Must specify species mass using 'speciesMass'");

    // Elementary charge
    if (tbl.hasNumber("speciesCharge"))
      speciesCharge = tbl.getNumber("speciesCharge");
    else
      throw Lucee::Except("BiasedSheath5DUpdater::readInput: Must specify species charge using 'speciesCharge'");
  }

  void
  BiasedSheath5DUpdater::initialize()
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

    // Fill out the node numbers on lower and upper surfaces in z
    lowerEdgeNodeNums = std::vector<int>(nodalBasis3d->getNumSurfLowerNodes(2));
    upperEdgeNodeNums = std::vector<int>(nodalBasis3d->getNumSurfUpperNodes(2));
    nodalBasis3d->getSurfLowerNodeNums(2, lowerEdgeNodeNums);
    nodalBasis3d->getSurfUpperNodeNums(2, upperEdgeNodeNums);
  }

  Lucee::UpdaterStatus
  BiasedSheath5DUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    const Lucee::Field<3, double>& phiIn = this->getInp<Lucee::Field<3, double> >(0);
    const Lucee::Field<5, double>& hamilDerivIn = this->getInp<Lucee::Field<5, double> >(1);
    // Output distribution function
    Lucee::Field<5, double>& distf = this->getOut<Lucee::Field<5, double> >(0);

    Lucee::Region<5, int> globalRgn = grid.getGlobalRegion();
    Lucee::Region<5, int> localRgn = grid.getLocalRegion();

    Lucee::ConstFieldPtr<double> phiInPtr = phiIn.createConstPtr();
    Lucee::ConstFieldPtr<double> hamilDerivInPtr = hamilDerivIn.createConstPtr();
    Lucee::FieldPtr<double> sknPtr = distf.createPtr(); // for skin-cell
    Lucee::FieldPtr<double> gstPtr = distf.createPtr(); // for ghost-cell

    unsigned nlocal = nodalBasis5d->getNumNodes();
    unsigned nlocal3d = nodalBasis3d->getNumNodes();
    double cellCentroid[5];
    // index of skin cell
    int idx[5];
    // index of ghost cell
    int gstIdx[5];

    // Apply sheath boundary conditions on lower surface
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
          phiIn.setPtr(phiInPtr, idx[0], idx[1], idx[2]);
          
          gstIdx[0] = ix;
          gstIdx[1] = iy;
          gstIdx[2] = idx[2]-1;
          // Loop over all configuration space nodes contained in this element
          for (int configNode = 0; configNode < lowerEdgeNodeNums.size(); configNode++)
          {
            int configNodeIndex = lowerEdgeNodeNums[configNode];

            // Do not do anything at this node if phi*charge >= 0
            if ( phiInPtr[configNodeIndex]*speciesCharge >= 0.0 )
              continue;

            // Cutoff velocity at this node
            double cutoffV = -std::sqrt(-2*phiInPtr[configNodeIndex]*speciesCharge/speciesMass);
            bool foundCutoffCell = false;

            // Need to loop over velocity space in this specific manner
            for (int ivSkin = localRgn.getLower(3), ivGhost = localRgn.getUpper(3)-1;
                ivSkin < localRgn.getUpper(3); ivSkin++, ivGhost--)
            {
              idx[3] = ivSkin;
              gstIdx[3] = ivGhost;

              // Check to see if there is no chance of matching the ion flux anymore
              idx[4] = localRgn.getLower(4);
              grid.setIndex(idx);
              grid.getCentroid(cellCentroid);
              // Skip everything if skin cell now has a positive parallel velocity
              if (cellCentroid[3] >= 0.0)
              {
                if (foundCutoffCell == false)
                {
                  std::cout << "(L) Unable to find a match to cutoff velocity" << std::endl;
                  std::cout << "charge = " << speciesCharge << std::endl;
                  std::cout << "cutoffV = " << cutoffV << std::endl;
                  std::cout << "cellLower = " << cellCentroid[3]-0.5*grid.getDx(3) << std::endl;
                  std::cout << "cellUpper = " << cellCentroid[3]+0.5*grid.getDx(3) << std::endl;
                }
                break;
              }
              // Check to see if cutoff v is above vmax
              if (cutoffV < -0.5*grid.getDx(3)*(globalRgn.getUpper(3)-globalRgn.getLower(3)))
              {
                // Reflect everything
                foundCutoffCell = true;
              }

              if (foundCutoffCell == false)
              {
                // Check if cutoffV falls in this cell
                if (cellCentroid[3] - 0.5*grid.getDx(3) <= cutoffV && 
                    cellCentroid[3] + 0.5*grid.getDx(3) >= cutoffV)
                {
                  foundCutoffCell = true;

                  // First compute total flux in this cell
                  double speciesFluxAtIv = 0.0;
                  // At every velocity space index, need to compute flux over all cells in mu
                  for (int iMu = localRgn.getLower(4); iMu < localRgn.getUpper(4); iMu++)
                  {
                    idx[4] = iMu;
                    distf.setPtr(sknPtr, idx);
                    hamilDerivIn.setPtr(hamilDerivInPtr, idx);
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
                      hamilReduced(nodeIndex) = hamilDerivInPtr[nodalStencil[nodeIndex] + configNodeIndex];
                    }
                    speciesFluxAtIv += scaleFactor*distfReduced.dot(momentMatrix*hamilReduced);
                  }

                  // Then compute outward flux only above cutoff velocity
                  double a = -0.5*grid.getDx(3);
                  double b = cutoffV - cellCentroid[3]; // Needs to be in local coordinates [-dv/2, dv/2]
                  // Integration weight scale for this modified integral
                  double weightScale = 0.5*(b-a)*0.5*grid.getDx(4);
                  double speciesFluxAboveVc = 0.0;
                  double refCoord[2];
                  std::vector<double> basisAtPoint(nodalStencil.size());
                  
                  for (int gaussIndex = 0; gaussIndex < gaussSurfCoords.rows(); gaussIndex++)
                  {
                    refCoord[0] = ( 0.5*(b-a)*gaussSurfCoords(gaussIndex,0) + 0.5*(a+b) )/(0.5*grid.getDx(3));
                    refCoord[1] = gaussSurfCoords(gaussIndex,1);
                    nodalBasis2d->evalBasis(refCoord, basisAtPoint);

                    for (int iMu = localRgn.getLower(4); iMu < localRgn.getUpper(4); iMu++)
                    {
                      idx[4] = iMu;
                      distf.setPtr(sknPtr, idx);
                      hamilDerivIn.setPtr(hamilDerivInPtr, idx);
                      // Compute distribution function at quadrature point
                      double fAtPoint = 0.0;
                      double gradHAtPoint = 0.0;
                      for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                      {
                        fAtPoint += sknPtr[nodalStencil[nodeIndex] + configNodeIndex]*basisAtPoint[nodeIndex];
                        gradHAtPoint += hamilDerivInPtr[nodalStencil[nodeIndex] + configNodeIndex]*basisAtPoint[nodeIndex];
                      }

                      speciesFluxAboveVc += weightScale*gaussSurfWeights(gaussIndex)*
                        scaleFactor*gradHAtPoint*fAtPoint;
                    }
                  }

                  // excessFraction should be between 0.0 and 1.0
                  double excessFraction = (speciesFluxAtIv - speciesFluxAboveVc)/speciesFluxAtIv;

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
                  // Make sure nothing is reflected at this velocity
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
                // Already iterated past cutoff cell, so just reflect everything
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
    // Apply sheath boundary conditions on lower surface
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
          phiIn.setPtr(phiInPtr, idx[0], idx[1], idx[2]);
          
          gstIdx[0] = ix;
          gstIdx[1] = iy;
          gstIdx[2] = idx[2]+1;
          // Loop over all configuration space nodes contained in this element
          for (int configNode = 0; configNode < upperEdgeNodeNums.size(); configNode++)
          {
            int configNodeIndex = upperEdgeNodeNums[configNode];

            // Do not do anything at this node if phi*charge > 0
            if ( phiInPtr[configNodeIndex]*speciesCharge >= 0.0 )
              continue;

            // Cutoff velocity at this node
            double cutoffV = std::sqrt(-2*phiInPtr[configNodeIndex]*speciesCharge/speciesMass);
            bool foundCutoffCell = false;

            // Need to loop over velocity space in this specific manner
            for (int ivSkin = localRgn.getUpper(3)-1, ivGhost = localRgn.getLower(3);
                ivSkin >= localRgn.getLower(3); ivSkin--, ivGhost++)
            {
              idx[3] = ivSkin;
              gstIdx[3] = ivGhost;

              // Check to see if there is no chance of matching the ion flux anymore
              idx[4] = localRgn.getLower(4);
              grid.setIndex(idx);
              grid.getCentroid(cellCentroid);
              // Skip everything if skin cell now has a negative parallel velocity
              if (cellCentroid[3] < 0.0)
              {
                if (foundCutoffCell == false)
                {
                  std::cout << "(U) Unable to find a match to cutoff velocity" << std::endl;
                  std::cout << "charge = " << speciesCharge << std::endl;
                  std::cout << "cutoffV = " << cutoffV << std::endl;
                  std::cout << "cellLower = " << cellCentroid[3]-0.5*grid.getDx(3) << std::endl;
                  std::cout << "cellUpper = " << cellCentroid[3]+0.5*grid.getDx(3) << std::endl;
                }
                break;
              }
              // Check to see if cutoff v is above vmax
              if (cutoffV > 0.5*grid.getDx(3)*(globalRgn.getUpper(3)-globalRgn.getLower(3)))
              {
                // Reflect everything
                foundCutoffCell = true;
              }

              if (foundCutoffCell == false)
              {
                // Check if cutoffV falls in this cell
                if (cellCentroid[3] - 0.5*grid.getDx(3) <= cutoffV && 
                    cellCentroid[3] + 0.5*grid.getDx(3) >= cutoffV)
                {
                  foundCutoffCell = true;

                  // First compute total flux in this cell
                  double speciesFluxAtIv = 0.0;
                  // At every velocity space index, need to compute flux over all cells in mu
                  for (int iMu = localRgn.getLower(4); iMu < localRgn.getUpper(4); iMu++)
                  {
                    idx[4] = iMu;
                    distf.setPtr(sknPtr, idx);
                    hamilDerivIn.setPtr(hamilDerivInPtr, idx);
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
                      hamilReduced(nodeIndex) = hamilDerivInPtr[nodalStencil[nodeIndex] + configNodeIndex];
                    }
                    speciesFluxAtIv += scaleFactor*distfReduced.dot(momentMatrix*hamilReduced);
                  }

                  // Then compute outward flux only above cutoff velocity
                  double a = cutoffV - cellCentroid[3]; // Needs to be in local coordinates [-dv/2, dv/2]
                  double b = 0.5*grid.getDx(3);
                  // Integration weight scale for this modified integral
                  double weightScale = 0.5*(b-a)*0.5*grid.getDx(4);
                  double speciesFluxAboveVc = 0.0;
                  double refCoord[2];
                  std::vector<double> basisAtPoint(nodalStencil.size());
                  
                  for (int gaussIndex = 0; gaussIndex < gaussSurfCoords.rows(); gaussIndex++)
                  {
                    refCoord[0] = ( 0.5*(b-a)*gaussSurfCoords(gaussIndex,0) + 0.5*(a+b) )/(0.5*grid.getDx(3));
                    refCoord[1] = gaussSurfCoords(gaussIndex,1);
                    nodalBasis2d->evalBasis(refCoord, basisAtPoint);

                    for (int iMu = localRgn.getLower(4); iMu < localRgn.getUpper(4); iMu++)
                    {
                      idx[4] = iMu;
                      distf.setPtr(sknPtr, idx);
                      hamilDerivIn.setPtr(hamilDerivInPtr, idx);
                      // Compute distribution function at quadrature point
                      double fAtPoint = 0.0;
                      double gradHAtPoint = 0.0;
                      for (int nodeIndex = 0; nodeIndex < nodalStencil.size(); nodeIndex++)
                      {
                        fAtPoint += sknPtr[nodalStencil[nodeIndex] + configNodeIndex]*basisAtPoint[nodeIndex];
                        gradHAtPoint += hamilDerivInPtr[nodalStencil[nodeIndex] + configNodeIndex]*basisAtPoint[nodeIndex];
                      }

                      speciesFluxAboveVc += weightScale*gaussSurfWeights(gaussIndex)*
                        scaleFactor*gradHAtPoint*fAtPoint;
                    }
                  }

                  // excessFraction should be between 0.0 and 1.0
                  double excessFraction = (speciesFluxAtIv - speciesFluxAboveVc)/speciesFluxAtIv;

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
                  // Make sure nothing is reflected at this velocity
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
                // Already iterated past cutoff cell, so just reflect everything
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

    return Lucee::UpdaterStatus();
  }

  void
  BiasedSheath5DUpdater::declareTypes()
  {
    // Potential (eV)
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Hamiltonian of species
    this->appendInpVarType(typeid(Lucee::Field<5, double>));
    // Distribution function of species to be reflected
    this->appendOutVarType(typeid(Lucee::Field<5, double>));
  }

  void
  BiasedSheath5DUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
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
  BiasedSheath5DUpdater::sameConfigCoords(int srcIndex, int tarIndex, double dxMin,
    const Eigen::MatrixXd& nodeList)
  {
    for (int d = 0; d < 3; d++)
      if (std::fabs(nodeList(srcIndex,d)-nodeList(tarIndex,d)) > 1e-4*dxMin) 
        return false;
    return true;
  }
}
