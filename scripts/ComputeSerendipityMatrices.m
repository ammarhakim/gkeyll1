% ComputeSerendipityMatrices.m
% Eric Shi (10-19-2012)
% Computes the serendipity basis functions for 4, 8, and 12 node elements
% in 2-D. Then evaluates the M, K, and G_i matrices.
% Needs some modifications to work in 3D.
clc
clear

dim       = 2;
% nodeList  = {[-1,-1],[1,-1],[1,1],[-1,1]};
% degree    = 1; % (r) This is tied to nodeList!
% nodeList  = {[-1,-1],[0,-1],[1,-1],[1,0],[1,1],[0,1],[-1,1],[-1,0]};
% degree    = 2; % (r) This is tied to nodeList!
nodeList  = {[-1,-1],[-0.5,-1],[0.5,-1],[1,-1],[1,-0.5],[1,0.5],[1,1],...
                [0.5,1],[-0.5,1],[-1,1],[-1,0.5],[-1,-0.5]};
degree    = 3; % (r) This is tied to nodeList!
basisList = {};

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

% Evaluate mass matrix
massMatrix = zeros(length(functionVector),length(functionVector));
for kIndex = 1:length(functionVector)
    for mIndex = 1:length(functionVector)
        f = functionVector(kIndex)*functionVector(mIndex);
        f = int(f, x, -1, 1);
        f = int(f, y, -1, 1);
        massMatrix(kIndex,mIndex) = f;
    end
end
% Put in readable format
massMatrix = sym(massMatrix)

% Evaluate K matrix
kMatrix = zeros(length(functionVector),length(functionVector));
jacobianMatrix = jacobian(functionVector); % First column is dx, second column is dy

for kIndex = 1:length(functionVector)
    for mIndex = 1:length(functionVector)
        f = jacobianMatrix(kIndex,1)*jacobianMatrix(mIndex,1) + ...
            jacobianMatrix(kIndex,2)*jacobianMatrix(mIndex,2);
        f = int(f, x, -1, 1);
        f = int(f, y, -1, 1);
        kMatrix(kIndex,mIndex) = f;
    end
end
% Put in readable format
kMatrix = sym(kMatrix)

% Evaluate G^{km}_x
gxMatrix = zeros(length(functionVector),length(functionVector));
for kIndex = 1:length(functionVector)
    for mIndex = 1:length(functionVector)
        f = jacobianMatrix(kIndex,1)*functionVector(mIndex);
        f = int(f, x, -1, 1);
        f = int(f, y, -1, 1);
        gxMatrix(kIndex,mIndex) = f;
    end
end
% Put in readable format
gxMatrix = sym(gxMatrix)

% Evaluate G^{km}_y
gyMatrix = zeros(length(functionVector),length(functionVector));
for kIndex = 1:length(functionVector)
    for mIndex = 1:length(functionVector)
        f = jacobianMatrix(kIndex,2)*functionVector(mIndex);
        f = int(f, x, -1, 1);
        f = int(f, y, -1, 1);
        gyMatrix(kIndex,mIndex) = f;
    end
end
% Put in readable format
gyMatrix = sym(gyMatrix)

% Output code in latex
% latex(functionVector)
% latex(massMatrix)
% latex(kMatrix)
% latex(gxMatrix)
% latex(gyMatrix)