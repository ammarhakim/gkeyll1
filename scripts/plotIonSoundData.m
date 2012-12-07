function growthRate = plotIonSoundData(numVPoints,tempRatioStr)

    fileDir = '/Users/dark1egion/Research/Gkeyll-Project/gkeyllall/ionsound/';
    fileBase = ['ionsound',num2str(numVPoints),'_',tempRatioStr,'_fieldEnergy_']

    fileNums = 0:20;

    timeVals = [];
    dataVals = [];

    for i = 1:length(fileNums)
        filename = [fileDir,fileBase,num2str(fileNums(i)),'.h5'];
        dataVec = h5read(filename,'/DataStruct/data');
        timeVec = h5read(filename,'/DataStruct/timeMesh');
        timeVals = [timeVals,timeVec];
        dataVals = [dataVals,dataVec];
    end

    dData = diff(dataVals); % Compute [x(2)-x(1),x(3)-x(2),..]

    % Find locations when slope changes from positive to negative
    extremeIndices = [];

    for i = 2:length(dData)
        if dData(i-1) >= 0 && dData(i) < 0
            extremeIndices = [extremeIndices,i];
        end
    end
    
    % Find first two adjacent extremes when energy decays, if it exists
    while dataVals(extremeIndices(1)) < dataVals(extremeIndices(2))
        extremeIndices = extremeIndices(2:end);
    end
    
    % Compute Least Squares Fit

    % If rsquare is not close enough to 1, remove data points from end until
    % desired rsquare is achieved

    [polyCoeffs,S] = polyfit(timeVals(extremeIndices),log(dataVals(extremeIndices)),1);
    rsquare = 1 - S.normr^2 / norm(log(dataVals(extremeIndices))-mean(log(dataVals(extremeIndices))))^2;

    while rsquare < 0.99 && length(extremeIndices) > 2
        extremeIndices = extremeIndices(1:end-1);
        [polyCoeffs,S] = polyfit(timeVals(extremeIndices),log(dataVals(extremeIndices)),1);
        rsquare = 1 - S.normr^2 / norm(log(dataVals(extremeIndices))-mean(log(dataVals(extremeIndices))))^2;
    end

    rsquare
    length(extremeIndices)

    extremeVals = dataVals(extremeIndices);

    lineX = timeVals;
    lineY = polyval(polyCoeffs,lineX);
    lineY = exp(lineY);

    hold on
    plot(timeVals,dataVals)
    plot(timeVals(extremeIndices),extremeVals,'.')
    plot(lineX,lineY)
    hold off
    ylim([min(extremeVals),max(extremeVals)]);

    growthRate = -polyCoeffs(1);

end