function plotSINDy(t,x,t2,x2)

figure
set(gcf,'Position',[75 75 450 750])

nP = 5;

subplot(nP,1,1)
plot(t,x(:,1),"Color",'b','LineWidth',1.2); hold on
plot(t2,x2(:,1),'k--')
ylabel('x')
legend({'Data','SINDy'})
set(gca,'XTickLabel',[]);
subplot(nP,1,2)
plot(t,x(:,2),"Color",'r','LineWidth',1.2); hold on
plot(t2,x2(:,2),'k--')
ylabel('y')
set(gca,'XTickLabel',[]);
subplot(nP,1,3)
plot(t,x(:,3),"Color",'g','LineWidth',1.2); hold on
plot(t2,x2(:,3),'k--')
ylabel('z')
xlabel('time')


subplot(nP,1,[4,5])
plot3(x(:,1),x(:,2),x(:,3),'b','LineWidth',0.5), hold on
plot3(x2(:,1),x2(:,2),x2(:,3),'k--','LineWidth',0.5), hold on
% plot3(x(k,1),x(k,2),x(k,3),'r.','LineWidth',2,'MarkerSize',10), hold off
axis([-40 40 -40 40 -10 50])
view(-140,20);
% axis off
% axis equal
xlabel('x')
ylabel('y')
zlabel('z')