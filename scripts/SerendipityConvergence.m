clear
close all

% Path to data files
serendipityPath = '/Users/dark1egion/Research/Gkeyll-Project/gkeyllall/gkeyll-tests/gkeyll-tests/tests/slvrs/serendipityConvergence/';
lagrangePath = '/Users/dark1egion/Research/Gkeyll-Project/gkeyllall/gkeyll-tests/gkeyll-tests/tests/slvrs/lagrangeConvergence/';

cellSizes = [4,8,16,32];

% Manual input data (could parse but...)
execTimesS = zeros(length(cellSizes),3);
execTimesS(:,1) = [0.546904 2.76353 18.4613 142.552];
execTimesS(:,2) = [1.42926 7.94137 55.3581 428.358];
execTimesS(:,3) = [2.60174 15.4653 110.327 856.953];

execTimesL = zeros(length(cellSizes),3);
execTimesL(:,1) = [0.546053 2.76703 18.4478 143.458];
execTimesL(:,2) = [2.00978 11.6285 81.9071 634.071];
execTimesL(:,3) = [7.04044 44.9869 340.633 2627.37];

errorMatrixS = ComputeMethodError(serendipityPath,cellSizes);
[errorMatrixL,lengthVector] = ComputeMethodError(lagrangePath,cellSizes);

%% Plot execution time vs error
set(0,'DefaultAxesFontSize',12);
% figure('OuterPosition',[50, 50, 1100, 800]);%,'DefaultAxesFontWeight','bold'
% Currently not plotting the third order Lagrange result on a 32x32x32 grid
% because results have too high an error--need to run at lower CFL number.
loglog(execTimesS(:,1),errorMatrixS(:,1),'g.-',...
    execTimesL(:,1),errorMatrixL(:,1),'g.--',...
    execTimesS(:,2),errorMatrixS(:,2),'r.-',...
    execTimesL(:,2),errorMatrixL(:,2),'b.-',...
    execTimesS(:,3),errorMatrixS(:,3),'r.--',...
    execTimesL(1:end-1,3),errorMatrixL(1:end-1,3),'b.--')
legend('Serendipity 1','Lagrange 1','Serendipity 2','Lagrange 2','Serendipity 3','Lagrange 3')
xlabel('Time (s)')
ylabel('Error')

outputExecTimes = [execTimesS execTimesL];
save serendipityTimeComparison.txt outputExecTimes -ASCII
dlmwrite('serendipityTimeComparison.txt',['# Rows are increasing grid size, Cols 1-3 Serendipity 1-3, Cols 4-6 Lagrange 1-3' 13 10 fileread('serendipityTimeComparison.txt')],'delimiter','');

outputErrorMatrices = [lengthVector execTimesS execTimesL];
save serendipityError.txt outputErrorMatrices -ASCII
dlmwrite('serendipityError.txt',['# Col 1 is deltaX, Cols 2-4 Serendipity 1-3 Errors, Cols 5-7 Lagrange 1-3 Errors' 13 10 fileread('serendipityError.txt')],'delimiter','');
