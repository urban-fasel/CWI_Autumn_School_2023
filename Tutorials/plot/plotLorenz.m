function plotLorenz(t,x)

figure
% set(gcf,'Position',[75 75 900 450])
set(gcf,'Position',[75 75 450 750])

nP = 5;

subplot(nP,1,1)
plot(t,x(:,1),"Color",'b')
ylabel('x')
set(gca,'XTickLabel',[]);
subplot(nP,1,2)
plot(t,x(:,2),"Color",'r')
ylabel('y')
set(gca,'XTickLabel',[]);
subplot(nP,1,3)
plot(t,x(:,3),"Color",'g')
ylabel('z')
xlabel('time')


subplot(nP,1,[4,5])
plot3(x(:,1),x(:,2),x(:,3),'k-','LineWidth',1), hold on
% plot3(x(k,1),x(k,2),x(k,3),'r.','LineWidth',2,'MarkerSize',10), hold off
axis([-40 40 -40 40 -10 50])
view(-140,20);
% axis off
% axis equal
xlabel('x')
ylabel('y')
zlabel('z')