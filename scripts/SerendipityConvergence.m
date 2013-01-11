clear
close all

% Path to data files
serendipityPath = '/Users/dark1egion/Research/Gkeyll-Project/gkeyllall/gkeyll-tests/gkeyll-tests/tests/slvrs/serendipityConvergence/';
lagrangePath = '/Users/dark1egion/Research/Gkeyll-Project/gkeyllall/gkeyll-tests/gkeyll-tests/tests/slvrs/lagrangeConvergence/';

cellSizes = [4,8,16,32,48];

% Manual input data (could parse but...)
execTimesS = zeros(length(cellSizes),4);
execTimesS(:,1) = [0.037081 0.360195 4.806 73.2142 387.605];
execTimesS(:,2) = [0.093965 1.0351 14.4054 219.66 1128.78];
execTimesS(:,3) = [0.252053 2.98708 42.3902 653.952 3335.11];
execTimesS(:,4) = [0.650139 8.40709 123.603 1908.57 9550.21];

execTimesL = zeros(length(cellSizes),4);
execTimesL(:,1) = [0.03828 0.360748 4.80756 73.1977 383.743];
execTimesL(:,2) = [0.254555 3.01086 43.0352 663.005 3324.48];
execTimesL(:,3) = [1.7516 22.9395 342.476 5352.65 0];
execTimesL(:,4) = [5.03932 69.3499 1059.42 0 0];

errorMatrixS = ComputeMethodError(serendipityPath,cellSizes,execTimesS);
[errorMatrixL,lengthVector] = ComputeMethodError(lagrangePath,cellSizes,execTimesL);

%% Plot execution time vs error
% set(0,'DefaultAxesFontSize',12,'DefaultAxesFontWeight','bold');
set(0,'DefaultAxesFontSize',24);
% figure('OuterPosition',[50, 50, 1100, 800]);%,'DefaultAxesFontWeight','bold'
figure('OuterPosition',[50, 50, 1000, 800]);
loglog(execTimesS(:,1),errorMatrixS(:,1),'g.-',...
    execTimesL(:,1),errorMatrixL(:,1),'g.--',...
    execTimesS(:,2),errorMatrixS(:,2),'r.-',...
    execTimesL(:,2),errorMatrixL(:,2),'r.--',...
    execTimesS(:,3),errorMatrixS(:,3),'b.-',...
    execTimesL(:,3),errorMatrixL(:,3),'b.--',...
    execTimesS(:,4),errorMatrixS(:,4),'m.-',...
    execTimesL(1:3,4),errorMatrixL(1:3,4),'m.--','MarkerSize',14)%,'LineWidth',1.5,'MarkerSize',10
legend('Serendipity 1','Lagrange 1','Serendipity 2','Lagrange 2','Serendipity 3','Lagrange 3','Serendipity 4','Lagrange 4')
xlabel('Time (s)')
ylabel('Error')
%axis tight
XL = xlim;
YL = ylim;
set(gcf,'PaperPositionMode','auto')
print -depsc2 serendipityConvergence.eps

% loglog(execTimesS(:,1),errorMatrixS(:,1),'g.-',...
%     execTimesL(:,1),errorMatrixL(:,1),'g.--','MarkerSize',14)%,'LineWidth',1.5,'MarkerSize',10
% legend('Serendipity 1','Lagrange 1')
% xlabel('Time (s)')
% ylabel('Error')
% xlim(XL);
% ylim(YL);
% set(gcf,'PaperPositionMode','auto')
% print -depsc2 serendipityConvergence1.eps
% 
% loglog(execTimesS(:,1),errorMatrixS(:,1),'g.-',...
%     execTimesL(:,1),errorMatrixL(:,1),'g.--',...
%     execTimesS(:,2),errorMatrixS(:,2),'r.-',...
%     execTimesL(:,2),errorMatrixL(:,2),'r.--','MarkerSize',14)%,'LineWidth',1.5,'MarkerSize',10
% legend('Serendipity 1','Lagrange 1','Serendipity 2','Lagrange 2')
% xlabel('Time (s)')
% ylabel('Error')
% xlim(XL);
% ylim(YL);
% set(gcf,'PaperPositionMode','auto')
% print -depsc2 serendipityConvergence2.eps
% 
% loglog(execTimesS(:,1),errorMatrixS(:,1),'g.-',...
%     execTimesL(:,1),errorMatrixL(:,1),'g.--',...
%     execTimesS(:,2),errorMatrixS(:,2),'r.-',...
%     execTimesL(:,2),errorMatrixL(:,2),'r.--',...
%     execTimesS(:,3),errorMatrixS(:,3),'b.-',...
%     execTimesL(:,3),errorMatrixL(:,3),'b.--','MarkerSize',14)%,'LineWidth',1.5,'MarkerSize',10
% legend('Serendipity 1','Lagrange 1','Serendipity 2','Lagrange 2','Serendipity 3','Lagrange 3')
% xlabel('Time (s)')
% ylabel('Error')
% xlim(XL);
% ylim(YL);
% set(gcf,'PaperPositionMode','auto')
% print -depsc2 serendipityConvergence3.eps

figure('OuterPosition',[50, 50, 1000, 800]);
loglog(lengthVector,errorMatrixS(:,1),'g.-',...
    lengthVector,errorMatrixL(:,1),'g.--',...
    lengthVector,errorMatrixS(:,2),'r.-',...
    lengthVector,errorMatrixL(:,2),'r.--',...
    lengthVector,errorMatrixS(:,3),'b.-',...
    lengthVector,errorMatrixL(:,3),'b.--',...
    lengthVector,errorMatrixS(:,4),'m.-',...
    lengthVector(1:3),errorMatrixL(1:3,4),'m.--','MarkerSize',14)
legend('Serendipity 1','Lagrange 1','Serendipity 2','Lagrange 2','Serendipity 3','Lagrange 3','Serendipity 4','Lagrange 4','Location','SouthEast')
xlabel('Element Size')
ylabel('Error')
axis tight
set(gcf,'PaperPositionMode','auto')
print -depsc2 serendipityError.eps