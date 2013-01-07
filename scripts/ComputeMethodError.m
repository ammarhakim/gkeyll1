function [errorMatrix,lengthVector,slopeVector] = ComputeMethodError(pathBase,cellSizes,execTimes)
% Compute the error of the "3d pulse" lua script output by summing up the
% absolute value of the differences between the initial and final solutions

    fileBase = {'order1_','order2_','order3_','order4_'};
    fileSuffix = '_q_';

    fileNums = [0,4];

    % Cols correspond to method order
    errorMatrix = zeros(length(cellSizes),length(fileBase));
%     slopeVector = zeros(length(fileBase),1);
    lengthVector = zeros(length(cellSizes),1);
    
    % Loop over the results from each order of a method
    for orderIndex = 1:length(fileBase)
        % Loop over the results of the simulation done on different cell
        % counts
        for sizeIndex = 1:length(cellSizes)
            if execTimes(sizeIndex,orderIndex) == 0
                continue
            end
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
            diffSum = diffSum/(size(dataStart,4)*size(dataStart,3)*size(dataStart,2)*cellSizes(orderIndex));

            errorMatrix(sizeIndex,orderIndex) = diffSum;
            lengthVector(sizeIndex) = globalDLengths(1);
        end

        % Determine slope of fit
        xPlot = [];
        yPlot = [];
        for itemIndex = 1:length(execTimes(orderIndex,:))
            if execTimes(orderIndex,itemIndex) == 0
                continue
            else
                xPlot = [xPlot,log(lengthVector(itemIndex))];
                yPlot = [yPlot,log(errorMatrix(itemIndex,orderIndex))];
            end
        end
        
        [polyCoeffs,S] = polyfit(xPlot,yPlot,1);
        slopeVector(orderIndex) = polyCoeffs(1);
        orderIndex
        rsquare = 1 - S.normr^2 / norm(yPlot-mean(yPlot))^2
        yFit = polyval(polyCoeffs,xPlot);
        
        plot(xPlot,yPlot,'.-',xPlot,yFit)
        w = waitforbuttonpress;
        close
         
    %     
    end
    slopeVector
end

