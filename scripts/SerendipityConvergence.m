clear
close all

% Path to data files
serendipityPath = '/Users/dark1egion/Research/Gkeyll-Project/gkeyllall/gkeyll-tests/gkeyll-tests/tests/slvrs/serendipityConvergence/';
lagrangePath = '/Users/dark1egion/Research/Gkeyll-Project/gkeyllall/gkeyll-tests/gkeyll-tests/tests/slvrs/lagrangeConvergence/';

cellSizes = [4,8,16,32];

% Manual input data (could parse but...)
execTimesS = zeros(length(cellSizes),3);
execTimesS(:,1) = [0.013989 0.123331 1.4644 22.591];
execTimesS(:,2) = [0.083322 0.931234 12.9819 199.278];
execTimesS(:,3) = [0.250375 2.97663 42.3626 654.312];

execTimesL = zeros(length(cellSizes),3);
execTimesL(:,1) = [0.016877 0.124448 1.47013 22.5999];
execTimesL(:,2) = [0.115783 1.37645 19.5869 301.046];
execTimesL(:,3) = [0.678743 8.69726 129.266 2016.38];

errorMatrixS = ComputeMethodError(serendipityPath,cellSizes);
[errorMatrixL,lengthVector] = ComputeMethodError(lagrangePath,cellSizes);

%% Plot execution time vs error
set(0,'DefaultAxesFontSize',12);
% figure('OuterPosition',[50, 50, 1100, 800]);%,'DefaultAxesFontWeight','bold'
loglog(execTimesS(:,1),errorMatrixS(:,1),'g.-',...
    execTimesL(:,1),errorMatrixL(:,1),'g.--',...
    execTimesS(:,2),errorMatrixS(:,2),'r.-',...
    execTimesL(:,2),errorMatrixL(:,2),'r.--',...
    execTimesS(:,3),errorMatrixS(:,3),'b.-',...
    execTimesL(1:end,3),errorMatrixL(1:end,3),'b.--')
legend('Serendipity 1','Lagrange 1','Serendipity 2','Lagrange 2','Serendipity 3','Lagrange 3')
xlabel('Time (s)')
ylabel('Error')
axis tight

% outputExecTimes = [execTimesS execTimesL];
% save serendipityTimeComparison.txt outputExecTimes -ASCII
% dlmwrite('serendipityTimeComparison.txt',['# Rows are increasing grid size, Cols 1-3 Serendipity 1-3, Cols 4-6 Lagrange 1-3' 13 10 fileread('serendipityTimeComparison.txt')],'delimiter','');
% 
% outputErrorMatrices = [lengthVector execTimesS execTimesL];
% save serendipityError.txt outputErrorMatrices -ASCII
% dlmwrite('serendipityError.txt',['# Col 1 is deltaX, Cols 2-4 Serendipity 1-3 Errors, Cols 5-7 Lagrange 1-3 Errors' 13 10 fileread('serendipityError.txt')],'delimiter','');
