function interpMatrix = ComputeSerendipityMatricesGeneral(dim, order)
% ComputeSerendipityMatrices3D.m
% James Juno (08-15-2015)
% Computes the serendipity basis functions for polynomial order up to 4
% in up to 3D. Then evaluates the interpoltion matrix.
%
% Inputs:   dim = dimension of the space you want to project to
%           order = polynomial order of the basis set
%
% Outputs:  interpMatrix = Interpolation matrix to be used to project
%           solution vectors on a finer grid

format long

nodeBase = linspace(-1, 1, order+1); %Base from which edge nodes will assigned
%Only used in 1D and 2D since 3D and beyond the pattern becomes more
%nontrivial
nodeList = {}; %List for nodal values
interpList = {}; %List of where we will interpolate to
basisList = {}; %List which contains the exponents for the monomials in our
%basis polynomials

%There might be a better way to do this, but since we won't be solving
%systems larger than 6D, it doesn't strike me as a big deal to figure out

if dim == 1
    syms x f
    for i=1:length(nodeBase)
        nodeList = [nodeList, {[nodeBase(i)]}];
    end
    celldisp(nodeList)
    
    for xIndex = 0:order
        basisList{end+1} = [xIndex];
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
    
    % Create the symbolic functions of x using the computed coefficients
    % and put them in functionVector
    coeffMatrixInv = sym(coeffMatrixInv);
    functionVector = sym(zeros(size(coeffMatrixInv,2),1));
    for basisIndex = 1:size(coeffMatrixInv,2)
        f = 0;
        for monomialIndex = 1:length(basisList)
            f = f + coeffMatrixInv(monomialIndex,basisIndex)*x^basisList{monomialIndex}(1);
        end
        functionVector(basisIndex) = f;
    end
    
    %Each entry of functionVector contains a polynomials which evaluates to
    %1 at the corresponding node (index value of functionVector corresponds
    %to index value in nodeList function is 1 at) and 0 at all the other
    %nodes. Now since we have the polynomial forms of the basis functions
    %we can solve for the interpolation matrix
    
    %The interpolation matrix I_{ij} = phi_j(x_i) 
    %The interpolation matrix has a number of columns equal to the number
    %of basis functions for that order and dimension and a number of rows
    %equal to the number of points we are interpolating to. It is not clear
    %to me if there is an "optimum" number of interpolation points so I
    %have followed the convention already in Gkeyll
    
    %The plotting routines will take the matrix solved for here
    %(coefficients will be hard coded because higher orders and higher
    %dimensions actually take time to calculate, so it does not seem worth
    %it to calculate the interpolation matrix on the fly) and then perform
    %f_interp = I*f where f is the vector of coefficients solved for by
    %Gkeyll for some important quantity (distribution function, moments,
    %etc.) and f_interp is the projected quantity we will use for plotting
    %and data analysis.
    
    if order == 1
        interpList = {-0.5, 0.5};
        
        interpMatrix = zeros(length(interpList),length(basisList));
        
        for i=1:size(interpMatrix, 1)
            for j =1:size(interpMatrix, 2)
                interpMatrix(i,j) = subs(functionVector(j), interpList(i));
            end
        end
        
    elseif order == 2
        interpList = {-2/3, 0, 2/3};
        
        interpMatrix = zeros(length(interpList),length(basisList));
        
        for i=1:size(interpMatrix, 1)
            for j =1:size(interpMatrix, 2)
                interpMatrix(i,j) = subs(functionVector(j), interpList(i));
            end
        end
        
    elseif order == 3
        interpList = {-3/4, -1/4, 1/4, 3/4};
        
        interpMatrix = zeros(length(interpList),length(basisList));
        
        for i=1:size(interpMatrix, 1)
            for j =1:size(interpMatrix, 2)
                interpMatrix(i,j) = subs(functionVector(j), interpList(i));
            end
        end
        
    elseif order == 4
        interpList = {-4/5, -2/5, 0, 2/5, 4/5};
        
        interpMatrix = zeros(length(interpList),length(basisList));
        
        for i=1:size(interpMatrix, 1)
            for j =1:size(interpMatrix, 2)
                interpMatrix(i,j) = subs(functionVector(j), interpList(i));
            end
        end
        
    end
    
elseif dim == 2
    syms x y f
    for i=1:length(nodeBase)
        nodeList = [nodeList, {[nodeBase(i), -1]}];
    end
    
    for i=2:length(nodeBase)-1
        nodeList = [nodeList, {[-1, nodeBase(i)]}];
        
        if order==4 & mod(i, (length(nodeBase)+1)/2)==0
            nodeList = [nodeList, {[nodeBase(i), nodeBase(i)]}];
        end
        
        nodeList = [nodeList, {[1, nodeBase(i)]}];
    end
    
    for i=1:length(nodeBase)
        nodeList = [nodeList, {[nodeBase(i), 1]}];
    end
    celldisp(nodeList)
    
    % Create basis monomials
    for yIndex = 0:order
        for xIndex = 0:order
            superDeg = (xIndex>1)*xIndex + (yIndex>1)*yIndex;
            if superDeg < order+1
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
    functionVector = sym(zeros(size(coeffMatrixInv,2),1));
    for basisIndex = 1:size(coeffMatrixInv,2)
        f = 0;
        for monomialIndex = 1:length(basisList)
            f = f + coeffMatrixInv(monomialIndex,basisIndex)*x^basisList{monomialIndex}(1)*y^basisList{monomialIndex}(2);
        end
        functionVector(basisIndex) = f;
    end
    
    %Each entry of functionVector contains a polynomials which evaluates to
    %1 at the corresponding node (index value of functionVector corresponds
    %to index value in nodeList function is 1 at) and 0 at all the other
    %nodes. Now since we have the polynomial forms of the basis functions
    %we can solve for the interpolation matrix
    
    %The interpolation matrix I_{ij} = phi_j(x_i) 
    %The interpolation matrix has a number of columns equal to the number
    %of basis functions for that order and dimension and a number of rows
    %equal to the number of points we are interpolating to. It is not clear
    %to me if there is an "optimum" number of interpolation points so I
    %have followed the convention already in Gkeyll
    
    %The plotting routines will take the matrix solved for here
    %(coefficients will be hard coded because higher orders and higher
    %dimensions actually take time to calculate, so it does not seem worth
    %it to calculate the interpolation matrix on the fly) and then perform
    %f_interp = I*f where f is the vector of coefficients solved for by
    %Gkeyll for some important quantity (distribution function, moments,
    %etc.) and f_interp is the projected quantity we will use for plotting
    %and data analysis.
    
    if order == 1
        
        order1Range = linspace(-0.5,0.5,2);
        for m=1:length(order1Range)
            for n=1:length(order1Range)
                interpList = [interpList, {[order1Range(n), order1Range(m)]}];
            end
        end

        interpMatrix = zeros(length(interpList),length(basisList));

        for i=1:size(interpMatrix, 1)
            for j =1:size(interpMatrix, 2)
                interpMatrix(i,j) = subs(functionVector(j),{x,y}, {interpList{i}(1),interpList{i}(2)});
            end
        end

    elseif order == 2
        
        order2Range = linspace(-2/3,2/3,3);
        for m=1:length(order2Range)
            for n=1:length(order2Range)
                interpList = [interpList, {[order2Range(n), order2Range(m)]}];
            end
        end
        
        interpMatrix = zeros(length(interpList), length(basisList));

        for i=1:size(interpMatrix, 1)
            for j =1:size(interpMatrix, 2)
                interpMatrix(i,j) = subs(functionVector(j),{x,y}, {interpList{i}(1),interpList{i}(2)});
            end
        end
        
    elseif order == 3
        
        order3Range = linspace(-3/4,3/4,4);
        for m=1:length(order3Range)
            for n=1:length(order3Range)
                interpList = [interpList, {[order3Range(n), order3Range(m)]}];
            end
        end
        
        interpMatrix = zeros(length(interpList), length(basisList));
        
        for i=1:size(interpMatrix, 1)
            for j =1:size(interpMatrix, 2)
                interpMatrix(i,j) = subs(functionVector(j),{x,y}, {interpList{i}(1),interpList{i}(2)});
            end
        end
        
    elseif order == 4
        
        order4Range = linspace(-4/5,4/5,5);
        for m=1:length(order4Range)
            for n=1:length(order4Range)
                interpList = [interpList, {[order4Range(n), order4Range(m)]}];
            end
        end
        celldisp(interpList)
        
        interpMatrix = zeros(length(interpList), length(basisList));
        
        for i=1:size(interpMatrix, 1)
            for j =1:size(interpMatrix, 2)
                interpMatrix(i,j) = subs(functionVector(j),{x,y}, {interpList{i}(1),interpList{i}(2)});
            end
        end
        
    end
    
elseif dim == 3
    syms x y z f
    
    if order == 1
        nodeList = {[-1,-1,-1], [1,-1,-1],...
                     [-1,1,-1],[1,1,-1],...
                     [-1,-1,1],[1,-1,1],...
                     [-1,1,1],[1,1,1]};
    elseif order == 2
        nodeList = {[-1,-1,-1],...
                    [0,-1,-1],...
                    [1,-1,-1],...
                    [-1,0,-1],...
                    [1,0,-1],...
                    [-1,1,-1],...
                    [0,1,-1],...
                    [1,1,-1],...
                    [-1,-1,0],...
                    [1,-1,0],...
                    [-1,1,0],...
                    [1,1,0],...
                    [-1,-1,1],...
                    [0,-1,1],...
                    [1,-1,1],...
                    [-1,0,1],...
                    [1,0,1],...
                    [-1,1,1],...
                    [0,1,1],...
                    [1,1,1]};
    elseif order == 3
        nodeList = {[-1,-1,-1],...
                    [-1.0/3.0,-1,-1],...
                    [1.0/3.0,-1,-1],...
                    [1,-1,-1],...
                    [-1,-1.0/3.0,-1],...
                    [1,-1.0/3.0,-1],...
                    [-1,1.0/3.0,-1],...
                    [1,1.0/3.0,-1],...
                    [-1,1,-1],...
                    [-1.0/3.0,1,-1],...
                    [1.0/3.0,1,-1],...
                    [1,1,-1],...
                    [-1,-1,-1.0/3.0],...
                    [1,-1,-1.0/3.0],...
                    [-1,1,-1.0/3.0],...
                    [1,1,-1.0/3.0],...
                    [-1,-1,1.0/3.0],...
                    [1,-1,1.0/3.0],...
                    [-1,1,1.0/3.0],...
                    [1,1,1.0/3.0],...
                    [-1,-1,1],...
                    [-1.0/3.0,-1,1],...
                    [1.0/3.0,-1,1],...
                    [1,-1,1],...
                    [-1,-1.0/3.0,1],...
                    [1,-1.0/3.0,1],...
                    [-1,1.0/3.0,1],...
                    [1,1.0/3.0,1],...
                    [-1,1,1],...
                    [-1.0/3.0,1,1],...
                    [1.0/3.0,1,1],...
                    [1,1,1]};
    elseif order == 4
        nodeList = {[-1,-1,-1],...
	                [-0.5,-1,-1],...
                    [0,-1,-1],...
                    [0.5,-1,-1],...
                    [1,-1,-1],...
                    [-1,-0.5,-1],...
                    [1,-0.5,-1],...
                    [-1,0,-1],...
	                [0,0,-1],...
                    [1,0,-1],...
                    [-1,0.5,-1],...
                    [1,0.5,-1],...
                    [-1,1,-1],...
	                [-0.5,1,-1],...
                    [0,1,-1],...
                    [0.5,1,-1],...
                    [1,1,-1],...
                    [-1,-1,-0.5],...
                    [1,-1,-0.5],...
                    [-1,1,-0.5],...
                    [1,1,-0.5],...
                    [-1,-1,0],...
	                [0,-1,0],...
                    [1,-1,0],...
                    [-1,0,0],...
                    [1,0,0],...
                    [-1,1,0],...
	                [0,1,0],...
                    [1,1,0],...
                    [-1,-1,0.5],...
                    [1,-1,0.5],...
                    [-1,1,0.5],...
                    [1,1,0.5],...
                    [-1,-1,1],...
	                [-0.5,-1,1],...
                    [0,-1,1],...
                    [0.5,-1,1],...
                    [1,-1,1],...
                    [-1,-0.5,1],...
                    [1,-0.5,1],...
                    [-1,0,1],...
	                [0,0,1],...
                    [1,0,1],...
                    [-1,0.5,1],...
                    [1,0.5,1],...
                    [-1,1,1],...
	                [-0.5,1,1],...
                    [0,1,1],...
                    [0.5,1,1],...
                    [1,1,1]};
    end
    
    celldisp(nodeList)
    
    % Create basis monomials
    for zIndex = 0:order
        for yIndex = 0:order
            for xIndex = 0:order
                superDeg = (xIndex>1)*xIndex + (yIndex>1)*yIndex+ (zIndex>1)*zIndex;
                if superDeg < order+1
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
    functionVector = sym(zeros(size(coeffMatrixInv,2),1));
    for basisIndex = 1:size(coeffMatrixInv,2)
        f = 0;
        for monomialIndex = 1:length(basisList)
            f = f + coeffMatrixInv(monomialIndex,basisIndex)*x^basisList{monomialIndex}(1)*y^basisList{monomialIndex}(2)*z^basisList{monomialIndex}(3);
        end
        functionVector(basisIndex) = f;
    end
    
    %Each entry of functionVector contains a polynomials which evaluates to
    %1 at the corresponding node (index value of functionVector corresponds
    %to index value in nodeList function is 1 at) and 0 at all the other
    %nodes. Now since we have the polynomial forms of the basis functions
    %we can solve for the interpolation matrix
    
    %The interpolation matrix I_{ij} = phi_j(x_i) 
    %The interpolation matrix has a number of columns equal to the number
    %of basis functions for that order and dimension and a number of rows
    %equal to the number of points we are interpolating to. It is not clear
    %to me if there is an "optimum" number of interpolation points so I
    %have followed the convention already in Gkeyll
    
    %The plotting routines will take the matrix solved for here
    %(coefficients will be hard coded because higher orders and higher
    %dimensions actually take time to calculate, so it does not seem worth
    %it to calculate the interpolation matrix on the fly) and then perform
    %f_interp = I*f where f is the vector of coefficients solved for by
    %Gkeyll for some important quantity (distribution function, moments,
    %etc.) and f_interp is the projected quantity we will use for plotting
    %and data analysis.
    
    if order == 1
        
        order1Range = linspace(-0.5,0.5,2);
        for l=1:length(order1Range)
            for m=1:length(order1Range)
                for n=1:length(order1Range)
                    interpList = [interpList, {[order1Range(n), order1Range(m), order1Range(l)]}];
                end
            end
        end
        
        interpMatrix = zeros(length(interpList),length(basisList));

        for i=1:size(interpMatrix, 1)
            for j =1:size(interpMatrix, 2)
                interpMatrix(i,j) = subs(functionVector(j),{x,y,z}, {interpList{i}(1),interpList{i}(2),interpList{i}(3)});
            end
        end

    elseif order == 2

        order2Range = linspace(-2/3,2/3,3);
        for l=1:length(order2Range)
            for m=1:length(order2Range)
                for n=1:length(order2Range)
                    interpList = [interpList, {[order2Range(n), order2Range(m), order2Range(l)]}];
                end
            end
        end
        
        interpMatrix = zeros(length(interpList), length(basisList));

        for i=1:size(interpMatrix, 1)
            for j =1:size(interpMatrix, 2)
                interpMatrix(i,j) = subs(functionVector(j),{x,y,z}, {interpList{i}(1),interpList{i}(2),interpList{i}(3)});
            end
        end
        
    elseif order == 3
        order3Range = linspace(-3/4,3/4,4);
        for l=1:length(order3Range)
            for m=1:length(order3Range)
                for n=1:length(order3Range)
                    interpList = [interpList, {[order3Range(n), order3Range(m), order3Range(l)]}];
                end
            end
        end
        
        interpMatrix = zeros(length(interpList), length(basisList));
        
        for i=1:size(interpMatrix, 1)
            for j =1:size(interpMatrix, 2)
                interpMatrix(i,j) = subs(functionVector(j),{x,y,z}, {interpList{i}(1),interpList{i}(2),interpList{i}(3)});
            end
        end
        
    elseif order == 4
        
        order4Range = linspace(-4/5,4/5,5);
        for l=1:length(order4Range)
            for m=1:length(order4Range)
                for n=1:length(order4Range)
                    interpList = [interpList, {[order4Range(n), order4Range(m), order4Range(l)]}];
                end
            end
        end
        
        interpMatrix = zeros(length(interpList), length(basisList));
        
        for i=1:size(interpMatrix, 1)
            for j =1:size(interpMatrix, 2)
                interpMatrix(i,j) = subs(functionVector(j),{x,y,z}, {interpList{i}(1),interpList{i}(2),interpList{i}(3)});
            end
        end
        
    end
    
    
% elseif dim == 4
%     syms x y z w f
% elseif dim == 5
%     syms x y z w v f
% elseif dim == 6
%     syms x y z w v u f
end

end