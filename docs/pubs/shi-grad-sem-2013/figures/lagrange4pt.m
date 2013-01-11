x = linspace(-1,1,1000);

xNodes = [-1,-1/3,1/3,1];
yNodes = zeros(1,length(xNodes));

l1 = (x-xNodes(2))./(xNodes(1)-xNodes(2)).*(x-xNodes(3))./(xNodes(1)-xNodes(3))...
        .*(x-xNodes(4))./(xNodes(1)-xNodes(4));
l2 = (x-xNodes(1))./(xNodes(2)-xNodes(1)).*(x-xNodes(3))./(xNodes(2)-xNodes(3))...
        .*(x-xNodes(4))./(xNodes(2)-xNodes(4));
l3 = (x-xNodes(1))./(xNodes(3)-xNodes(1)).*(x-xNodes(2))./(xNodes(3)-xNodes(2))...
        .*(x-xNodes(4))./(xNodes(3)-xNodes(4));
l4 = (x-xNodes(1))./(xNodes(4)-xNodes(1)).*(x-xNodes(2))./(xNodes(4)-xNodes(2))...
    .*(x-xNodes(3))./(xNodes(4)-xNodes(3));

set(0,'DefaultAxesFontSize',34,'DefaultAxesFontWeight','bold');
figure('OuterPosition',[50, 50, 1000, 800]);
hold on
%line([-1 1],[0 0],'Color',[0 0 0],'LineWidth',3)
plot(x,l1,x,l2,x,l3,x,l4,'LineWidth',3)
plot(xNodes,yNodes,'.-k','LineWidth',3,'MarkerSize',30)
hold off
axis tight

set(gcf,'PaperPositionMode','auto')
print -depsc2 lagrange4.eps