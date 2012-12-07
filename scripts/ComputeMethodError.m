function [errorMatrix,lengthVector] = ComputeMethodError(pathBase,cellSizes)
% Compute the error of the "3d pulse" lua script output by summing up the
% absolute value of the differences between the initial and final solutions

    fileBase = {'order1_','order2_','order3_'};
    fileSuffix = '_q_';

    fileNums = [0,4];
    
    nodesPerCell = [8,20,32];

    % Cols correspond to method order
    errorMatrix = zeros(length(cellSizes),length(fileBase));
%     slopeVector = zeros(length(fileBase),1);
    lengthVector = zeros(length(cellSizes),1);
    
    % Loop over the results from each order of a method
    for orderIndex = 1:length(fileBase)
        % Loop over the results of the simulation done on different cell
        % counts
        for sizeIndex = 1:length(cellSizes)
            filename = [pathBase,fileBase{orderIndex},num2str(cellSizes(sizeIndex)),fileSuffix,num2str(fileNums(1)),'.h5'];
            vsNumCells    = double(h5readatt(filename,'/StructGrid','vsNumCells'));
            vsUpperBounds = h5readatt(filename,'/StructGrid','vsUpperBounds');
            vsLowerBounds = h5readatt(filename,'/StructGrid','vsLowerBounds');
            globalSizes = double(vsUpperBounds - vsLowerBounds);
            globalDLengths = globalSizes./vsNumCells;
            dataStart = h5read(filename,'/StructGridField');

            filename = [pathBase,fileBase{orderIndex},num2str(cellSizes(sizeIndex)),fileSuffix,num2str(fileNums(2)),'.h5'];
            dataEnd   = h5read(filename,'/StructGridField');

            % Loop over 4-d array to sum up elements
            diffSum = 0.0;

            for xIndex = 1:size(dataStart,4)
                for yIndex = 1:size(dataStart,3)
                    for zIndex = 1:size(dataStart,2)
                        diffSum = diffSum + sum(abs(dataStart(:,zIndex,yIndex,xIndex) - dataEnd(:,zIndex,yIndex,xIndex)));
                    end
                end
            end

            % Normalize by dividing by Nx*Ny*Nz*nodes
            diffSum = diffSum/(size(dataStart,4)*size(dataStart,3)*size(dataStart,2)*nodesPerCell(orderIndex));

            errorMatrix(sizeIndex,orderIndex) = diffSum;
            lengthVector(sizeIndex) = globalDLengths(1);
        end

        % Determine slope of fit
%         xPlot = log(lengthVector);
%         yPlot = log(errorMatrix(:,orderIndex));
%         [polyCoeffs,S] = polyfit(xPlot,yPlot,1);
%         %yFit = polyval(polyCoeffs,xPlot);
%         slopeVector(orderIndex) = polyCoeffs(1);
    %     plot(xPlot,yPlot,'.-',xPlot,yFit)
    end

end

