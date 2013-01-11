x = linspace(-1,1,1000);

xNodes = [-1,0,1];
yNodes = zeros(1,length(xNodes));

lm = (x-xNodes(2))./(xNodes(1)-xNodes(2)).*(x-xNodes(3))./(xNodes(1)-xNodes(3));
lz = (x-xNodes(1))./(xNodes(2)-xNodes(1)).*(x-xNodes(3))./(xNodes(2)-xNodes(3));
lp = (x-xNodes(1))./(xNodes(3)-xNodes(1)).*(x-xNodes(2))./(xNodes(3)-xNodes(2));

set(0,'DefaultAxesFontSize',40,'DefaultAxesFontWeight','bold');
figure('OuterPosition',[50, 50, 1000, 800]);
hold on
%line([-1 1],[0 0],'Color',[0 0 0],'LineWidth',3)
plot(x,lm,x,lz,x,lp,'LineWidth',3)
plot(xNodes,yNodes,'.-k','LineWidth',3,'MarkerSize',30)
hold off
axis tight
ylim([-0.3156    1.0563])
set(gcf,'PaperPositionMode','auto')
print -depsc2 lagrange3.eps