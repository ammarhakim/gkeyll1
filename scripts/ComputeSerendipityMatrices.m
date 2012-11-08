% ComputeSerendipityMatrices.m
% Eric Shi (10-19-2012)
% Computes the serendipity basis functions for 4, 8, and 12 node elements
% in 2-D. Then evaluates the M, K, and G_i matrices.
% Needs some modifications to work in 3D.
clc
clear
close all

% 2-D Arbitrary Quadrilateral to map reference element
% Prefer first element be lower left corner and going around CCW
quadList  = {[-2,-1],[1.5,-2],[2,1],[-1,2]};

dim       = 2;
% nodeList  = {[-1,-1],[1,-1],[1,1],[-1,1]};
% degree    = 1; % (r) This is tied to nodeList!
nodeList  = {[-1,-1],[0,-1],[1,-1],[1,0],[1,1],[0,1],[-1,1],[-1,0]};
degree    = 2; % (r) This is tied to nodeList!
% nodeList  = {[-1,-1],[-1/3,-1],[1/3,-1],[1,-1],[1,-1/3],[1,1/3],[1,1],...
%                 [1/3,1],[-1/3,1],[-1,1],[-1,1/3],[-1,-1/3]};
% degree    = 3; % (r) This is tied to nodeList!
basisList = {};

% Plot the reference element with nodes
hold on
for i = 1:length(nodeList)
    plot(nodeList{i}(1),nodeList{i}(2),'b*')
end
% Plot the corner nodes of the 2-D arbitrary quadrilateral
for i = 1:length(quadList)
    plot(quadList{i}(1),quadList{i}(2),'r*')
end
axis square
hold off

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
for yIndex = 0:degree
    for xIndex = 0:degree
        superDeg = (xIndex>1)*xIndex + (yIndex>1)*yIndex;
        if superDeg < degree+1
            basisList{end+1} = [xIndex,yIndex,superDeg];
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

% Create the symbolic functions of x and y using the computed coefficients
% and put them in functionVector
coeffMatrixInv = sym(coeffMatrixInv);
syms x y f
functionVector = sym(zeros(size(coeffMatrixInv,2),1));
for basisIndex = 1:size(coeffMatrixInv,2)
    f = 0;
    for monomialIndex = 1:length(basisList)
        f = f + coeffMatrixInv(monomialIndex,basisIndex)*x^basisList{monomialIndex}(1)*y^basisList{monomialIndex}(2);
    end
    functionVector(basisIndex) = f;
end

% Make sure basis functions evaluate to 1 at a node and 0 elsewhere by
% evaluating the functions defined in functionVector at each node.
% testMatrix = zeros(length(nodeList),length(functionVector));
% for basisIndex = 1:length(functionVector)
%     for nodeIndex = 1:length(nodeList)
%         testMatrix(nodeIndex, basisIndex) = subs(functionVector(basisIndex),[x,y],[nodeList{nodeIndex}(1),nodeList{nodeIndex}(2)]);
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
jacobianDet   = det(jacobian(transformationMatrix)); % Jacobian of the transformation matrix
% Should name the following variable something else to reduce confusion
% with above variable
jacobianMatrix = jacobian(functionVector); % First column is dx, second column is dy

% Evaluate mass matrix
massMatrix = zeros(length(functionVector),length(functionVector));
for kIndex = 1:length(functionVector)
    for mIndex = 1:length(functionVector)
        f = functionVector(kIndex)*functionVector(mIndex)*jacobianDet;
        f = int(f, x, -1, 1);
        f = int(f, y, -1, 1);
        massMatrix(kIndex,mIndex) = f;
    end
end

% Evaluate K matrix
kMatrix = zeros(length(functionVector),length(functionVector));
for kIndex = 1:length(functionVector)
    for mIndex = 1:length(functionVector)
        f = jacobianMatrix(kIndex,1)*jacobianMatrix(mIndex,1) + ...
            jacobianMatrix(kIndex,2)*jacobianMatrix(mIndex,2);
        f = jacobianDet*f;
        f = int(f, x, -1, 1);
        f = int(f, y, -1, 1);
        kMatrix(kIndex,mIndex) = f;
    end
end

% Evaluate G^{km}_x
gxMatrix = zeros(length(functionVector),length(functionVector));
for kIndex = 1:length(functionVector)
    for mIndex = 1:length(functionVector)
        f = jacobianMatrix(kIndex,1)*functionVector(mIndex)*jacobianDet;
        f = int(f, x, -1, 1);
        f = int(f, y, -1, 1);
        gxMatrix(kIndex,mIndex) = f;
    end
end

% Evaluate G^{km}_y
gyMatrix = zeros(length(functionVector),length(functionVector));
for kIndex = 1:length(functionVector)
    for mIndex = 1:length(functionVector)
        f = jacobianMatrix(kIndex,2)*functionVector(mIndex)*jacobianDet;
        f = int(f, x, -1, 1);
        f = int(f, y, -1, 1);
        gyMatrix(kIndex,mIndex) = f;
    end
end

% Put in readable format
massMatrix = sym(massMatrix)
kMatrix = sym(kMatrix)
gxMatrix = sym(gxMatrix)
gyMatrix = sym(gyMatrix)

% Generate output code in latex
% latex(functionVector)
% latex(massMatrix)
% latex(kMatrix)
% latex(gxMatrix)
% latex(gyMatrix)