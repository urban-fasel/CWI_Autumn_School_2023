function animateLorenz(t,x)

T = 200; % plot every 100 steps
L = length(t);

figure
% set(gcf,'Position',[75 75 900 450])
set(gcf,'Position',[75 75 450 750])

nP1 = 5;%3
nP2 = 1;%2

for k=2:T:L
subplot(nP1,nP2,1)
plot(t(1:k,1),x(1:k,1),"Color",'b')
% plot(t,x(k,1),"Color",'b')
ylabel('x')
xlim([0 t(end)])
ylim([-40 40])
subplot(nP1,nP2,2)
plot(t(1:k,1),x(1:k,2),"Color",'r')
ylabel('y')
xlim([0 t(end)])
ylim([-40 40])
subplot(nP1,nP2,3)
plot(t(1:k,1),x(1:k,3),"Color",'g')
ylabel('z')
xlabel('time')
xlim([0 t(end)])
ylim([-10 50])


% subplot(nP1,nP2,[2,4,6])
subplot(nP1,nP2,[4,5])
plot3(x(1:k,1),x(1:k,2),x(1:k,3),'k-','LineWidth',1), hold on
plot3(x(k,1),x(k,2),x(k,3),'r.','LineWidth',2,'MarkerSize',10), hold off
axis([-40 40 -40 40 -10 50])
view(-140,20);
view(-140+230*k/L,24);
% axis off
% axis equal
xlabel('x')
ylabel('y')
zlabel('z')

drawnow

end