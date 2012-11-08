% ComputeSerendipityMatrices3D.m
% Eric Shi (11-03-2012)
% Computes the serendipity basis functions for 8, 20, and 32 node 'bricks'
% in 3-D. Then evaluates the M, K, and G_i matrices.
clc
clear
close all

% Corners of the 3-D Arbitrary six faced element
% quadList  = {[-2,-1,-1],[1.5,-2,-2],[2,1,0],[-1,2,0],...
%              [-2,-1,2],[1.5,-2,1],[2,1,3],[-1,2,1]};
quadList  = {[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],...
             [-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]};

dim       = 3;

% nodeList  = {[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],...
%              [-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]};
% degree    = 1; % (r) This is tied to nodeList!
% nodeList  = {[-1,-1,-1],[0,-1,-1],[1,-1,-1],[1,0,-1],[1,1,-1],[0,1,-1],[-1,1,-1],[-1,0,-1],...% End of Bottom Layer
%              [-1,-1,0],[1,-1,0],[1,1,0],[-1,1,0],...% End of Middle Layer
%              [-1,-1,1],[0,-1,1],[1,-1,1],[1,0,1],[1,1,1],[0,1,1],[-1,1,1],[-1,0,1]};% End of Top Layer
% degree    = 2; % (r) This is tied to nodeList!
nodeList  = {[-1,-1,-1],[-1/3,-1,-1],[1/3,-1,-1],[1,-1,-1],[1,-1/3,-1],[1,1/3,-1],[1,1,-1],...
                [1/3,1,-1],[-1/3,1,-1],[-1,1,-1],[-1,1/3,-1],[-1,-1/3,-1],... End of Bottom Layer
                [-1,-1,-1/3],[1,-1,-1/3],[1,1,-1/3],[-1,1,-1/3],... % End of 2nd Layer
                [-1,-1,1/3],[1,-1,1/3],[1,1,1/3],[-1,1,1/3],... % 3rd Layer
                [-1,-1,1],[-1/3,-1,1],[1/3,-1,1],[1,-1,1],[1,-1/3,1],[1,1/3,1],[1,1,1],...
                [1/3,1,1],[-1/3,1,1],[-1,1,1],[-1,1/3,1],[-1,-1/3,1]}; % End of Top Layer
degree    = 3; % (r) This is tied to nodeList!
basisList = {};

% % Plot the Nodes in 3D
% xPoints = zeros(length(nodeList),1);
% yPoints = zeros(length(nodeList),1);
% zPoints = zeros(length(nodeList),1);
% for i = 1:length(nodeList)
%     xPoints(i) = nodeList{i}(1);
%     yPoints(i) = nodeList{i}(2);
%     zPoints(i) = nodeList{i}(3);
% end
% xPointsArb = zeros(length(quadList),1);
% yPointsArb = zeros(length(quadList),1);
% zPointsArb = zeros(length(quadList),1);
% for i = 1:length(quadList)
%     xPointsArb(i) = quadList{i}(1);
%     yPointsArb(i) = quadList{i}(2);
%     zPointsArb(i) = quadList{i}(3);
% end
% % Plot the Nodes in 3D
% hold on
% scatter3(xPoints,yPoints,zPoints,'b')
% scatter3(xPointsArb,yPointsArb,zPointsArb,'r')
% hold off

% Figure out where the corner nodes are (arbitrary dimension)
% put these indices in the vector refCorners
distanceVector = zeros(length(nodeList),1); % Distance from origin squared
for nodeIndex = 1:length(nodeList)
    tempNode = nodeList{nodeIndex};
    distanceVector(nodeIndex) = sum(tempNode.^2);
end
maxDistance = max(distanceVector);
refCorners = zeros(length(quadList),1);
refIndex = 1;
for nodeIndex = 1:length(nodeList)
    if sum(nodeList{nodeIndex}.^2) == maxDistance
        refCorners(refIndex) = nodeIndex;
        refIndex = refIndex + 1;
    end
end

% Create basis monomials
for zIndex = 0:degree
    for yIndex = 0:degree
        for xIndex = 0:degree
            superDeg = (xIndex>1)*xIndex + (yIndex>1)*yIndex + (zIndex>1)*zIndex;
            if superDeg < degree+1
                basisList{end+1} = [xIndex,yIndex,zIndex,superDeg];
            end
        end
    end
end

coeffMatrix = zeros(length(nodeList), length(basisList));

for nodeIndex = 1:length(nodeList)
    nodePos = nodeList{nodeIndex};
    for basisIndex = 1:length(basisList)
        % Want to compute the monomial terms evaluated at a particular
        % coordinate
        coeffProduct = 1;
        basisTerm = basisList{basisIndex};
        for dimIndex = 1:dim
            coeffProduct = coeffProduct*nodePos(dimIndex)^basisTerm(dimIndex);
        end
        coeffMatrix(nodeIndex,basisIndex) = coeffProduct;
    end
end

coeffMatrixInv = inv(coeffMatrix);

% Create the symbolic functions of x, y, and z using the computed coefficients
% and put them in functionVector
coeffMatrixInv = sym(coeffMatrixInv);
syms x y z f
functionVector = sym(zeros(size(coeffMatrixInv,2),1));
for basisIndex = 1:size(coeffMatrixInv,2)
    f = 0;
    for monomialIndex = 1:length(basisList)
        f = f + coeffMatrixInv(monomialIndex,basisIndex)*x^basisList{monomialIndex}(1)*y^basisList{monomialIndex}(2)*z^basisList{monomialIndex}(3);
    end
    functionVector(basisIndex) = f;
end

% % Make sure basis functions evaluate to 1 at a node and 0 elsewhere by
% % evaluating the functions defined in functionVector at each node.
% testMatrix = zeros(length(nodeList),length(functionVector));
% for basisIndex = 1:length(functionVector)
%     for nodeIndex = 1:length(nodeList)
%         testMatrix(nodeIndex, basisIndex) = subs(functionVector(basisIndex),[x,y,z],[nodeList{nodeIndex}(1),nodeList{nodeIndex}(2),nodeList{nodeIndex}(3)]);
%     end
% end

% Create transformation x(eta,nu) and y(eta,nu) between reference nodes
% (nodeList) and quadList
transformationMatrix = sym(zeros(dim,1));
for nodeIndex = 1:length(quadList)
    for dimIndex = 1:dim
        transformationMatrix(dimIndex,1) = quadList{nodeIndex}(dimIndex)*functionVector(refCorners(nodeIndex)) + transformationMatrix(dimIndex,1);
    end
end
jacobianDet = det(jacobian(transformationMatrix));
% Should name the following variable something else to reduce confusion
% with above variable
jacobianMatrix = jacobian(functionVector); % First column is dx, second column is dy, third is dz

% % Evaluate mass matrix
% massMatrix = zeros(length(functionVector),length(functionVector));
% for kIndex = 1:length(functionVector)
%     for mIndex = 1:length(functionVector)
%         f = functionVector(kIndex)*functionVector(mIndex)*jacobianDet;
%         f = int(f, x, -1, 1);
%         f = int(f, y, -1, 1);
%         f = int(f, z, -1, 1);
%         massMatrix(kIndex,mIndex) = f;
%     end
% end

% % Evaluate K matrix
% kMatrix = zeros(length(functionVector),length(functionVector));
% for kIndex = 1:length(functionVector)
%     for mIndex = 1:length(functionVector)
%         f = jacobianMatrix(kIndex,1)*jacobianMatrix(mIndex,1) + ...
%             jacobianMatrix(kIndex,2)*jacobianMatrix(mIndex,2) + ...
%             jacobianMatrix(kIndex,3)*jacobianMatrix(mIndex,3);
%         f = jacobianDet*f;
%         f = int(f, x, -1, 1);
%         f = int(f, y, -1, 1);
%         f = int(f, z, -1, 1);
%         kMatrix(kIndex,mIndex) = f;
%     end
% end
% 
% % Evaluate G^{km}_x
% gxMatrix = zeros(length(functionVector),length(functionVector));
% for kIndex = 1:length(functionVector)
%     for mIndex = 1:length(functionVector)
%         f = jacobianMatrix(kIndex,1)*functionVector(mIndex)*jacobianDet;
%         f = int(f, x, -1, 1);
%         f = int(f, y, -1, 1);
%         f = int(f, z, -1, 1);
%         gxMatrix(kIndex,mIndex) = f;
%     end
% end
% 
% % Evaluate G^{km}_y
% gyMatrix = zeros(length(functionVector),length(functionVector));
% for kIndex = 1:length(functionVector)
%     for mIndex = 1:length(functionVector)
%         f = jacobianMatrix(kIndex,2)*functionVector(mIndex)*jacobianDet;
%         f = int(f, x, -1, 1);
%         f = int(f, y, -1, 1);
%         f = int(f, z, -1, 1);
%         gyMatrix(kIndex,mIndex) = f;
%     end
% end
% 
% % Evaluate G^{km}_z
% gzMatrix = zeros(length(functionVector),length(functionVector));
% for kIndex = 1:length(functionVector)
%     for mIndex = 1:length(functionVector)
%         f = jacobianMatrix(kIndex,3)*functionVector(mIndex)*jacobianDet;
%         f = int(f, x, -1, 1);
%         f = int(f, y, -1, 1);
%         f = int(f, z, -1, 1);
%         gzMatrix(kIndex,mIndex) = f;
%     end
% end

% % Put in readable format
% massMatrix = sym(massMatrix)
% kMatrix = sym(kMatrix)
% gxMatrix = sym(gxMatrix)
% gyMatrix = sym(gyMatrix)
% gzMatrix = sym(gzMatrix)

% Generate output code in latex
% latex(functionVector)
% latex(massMatrix)
% latex(kMatrix)
% latex(gxMatrix)
% latex(gyMatrix)