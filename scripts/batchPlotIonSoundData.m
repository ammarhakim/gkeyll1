% Goes through the various grid points and temp ratios and calculates
% the instability growth rates

% Add path to directory will data files
addpath('Users/dark1egion/Research/Gkeyll-Project/gkeyllall/ionsound/serendipityTests')

tRatios = {'0.1','0.3','0.5','0.75','1.0','1.5','2.0'};
% Data files for 8 have been created using a higher polyOrder
% 8,16,32,64,128
resolutionList = [8 16 32];
% First column is reserved for T_ratios
growthRateList = zeros(length(tRatios)+1,length(resolutionList)+1);

% Replace '.' with '_'
for i = 1:length(tRatios)
    % Store the t-ratio in the first column of the output structure
    growthRateList(i+1,1) = str2double(tRatios{i});
    tRatios{i} = strrep(tRatios{i},'.','_');
end

for j = 1:length(resolutionList)
    % Store the v-resolution in the first row of the appropriate column
    growthRateList(1,j+1) = resolutionList(j);
    for i = 1:length(tRatios)
        growthRateList(i+1,j+1) = plotIonSoundData(resolutionList(j),tRatios{i});
        w = waitforbuttonpress;
        close
    end
end

save simGrowthRates.txt growthRateList -ASCII