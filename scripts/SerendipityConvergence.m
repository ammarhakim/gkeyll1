clear
close all

% Path to data files
serendipityPath = '/Users/dark1egion/Research/Gkeyll-Project/gkeyllall/gkeyll-tests/gkeyll-tests/tests/slvrs/serendipityConvergence/';
lagrangePath = '/Users/dark1egion/Research/Gkeyll-Project/gkeyllall/gkeyll-tests/gkeyll-tests/tests/slvrs/lagrangeConvergence/';

cellSizes = [4,8,16,32];

% Manual input data (could parse but...)
execTimesS = zeros(length(cellSizes),3);
execTimesS(:,1) = [0.03081 0.360195 4.806 73.2142];
execTimesS(:,2) = [0.093965 1.0351 14.4054 219.66];
execTimesS(:,3) = [0.252053 2.98708 42.3902 653.952];

execTimesL = zeros(length(cellSizes),3);
execTimesL(:,1) = [0.03828 0.360748 4.80756 73.1977];
execTimesL(:,2) = [0.254555 3.01086 43.0352 663.005];
execTimesL(:,3) = [1.7516 22.9395 342.476 5352.65];

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
    execTimesL(:,3),errorMatrixL(:,3),'b.--')
legend('Serendipity 1','Lagrange 1','Serendipity 2','Lagrange 2','Serendipity 3','Lagrange 3')
xlabel('Time (s)')
ylabel('Error')
axis tight

figure
loglog(lengthVector,errorMatrixS(:,1),'g.-',...
    lengthVector,errorMatrixL(:,1),'g.--',...
    lengthVector,errorMatrixS(:,2),'r.-',...
    lengthVector,errorMatrixL(:,2),'r.--',...
    lengthVector,errorMatrixS(:,3),'b.-',...
    lengthVector,errorMatrixL(:,3),'b.--')
legend('Serendipity 1','Lagrange 1','Serendipity 2','Lagrange 2','Serendipity 3','Lagrange 3','Location','SouthEast')
xlabel('Element Size')
ylabel('Error')
axis tight

% outputExecTimes = [execTimesS execTimesL];
% save serendipityTimeComparison.txt outputExecTimes -ASCII
% dlmwrite('serendipityTimeComparison.txt',['# Rows are increasing grid size, Cols 1-3 Serendipity 1-3, Cols 4-6 Lagrange 1-3' 13 10 fileread('serendipityTimeComparison.txt')],'delimiter','');
% 
% outputErrorMatrices = [lengthVector execTimesS execTimesL];
% save serendipityError.txt outputErrorMatrices -ASCII
% dlmwrite('serendipityError.txt',['# Col 1 is deltaX, Cols 2-4 Serendipity 1-3 Errors, Cols 5-7 Lagrange 1-3 Errors' 13 10 fileread('serendipityError.txt')],'delimiter','');
