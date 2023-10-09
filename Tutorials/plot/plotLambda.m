function plotLambda(Theta,dx,lambda,n,m)

XiPath = [];
sumPath = [];

% phi = linspace(0,2*pi,m*n*2)';
phi = linspace(0,2*pi,m*n)';
phi = [phi;0];

dP = 0.2; % plot factor non active terms

for i = 1:length(lambda)
    Xi = sparsifyDynamics(Theta,dx,lambda(i),n); % identify model coefficients
    
    nonzero = Xi ~= 0;
    nonzero = double(nonzero(:));
    sumPath(i) = sum(nonzero);
    nonzero(nonzero==0)=dP;

    nn(i) = norm(Theta*Xi-dx);

    % nonzero = [nonzero,nonzero]';
    % nonzero = nonzero(:);
    nonzero = [nonzero;nonzero(1)];

    XiPath(:,i) = nonzero; 
end


figure
set(gcf,'Position',[0 0 400 200])

nat = 7; % number of active terms

% for i = 1:length(lambda)

% subplot(2,1,1)
semilogx(lambda,sumPath,'k'), hold on
semilogx([lambda(2) lambda(end)],[nat nat],'b--'), hold on
% semilogx(lambda(i),sumPath(i),'r*'), hold off
ylim([0 m*n])
xlabel('\lambda')
ylabel('number of active terms')

legend({'SINDy path', 'Lorenz system number of terms'})


% figure
% set(gcf,'Position',[0 0 400 200])
% 
% nat = 7; % number of active terms
% 
% % for i = 1:length(lambda)
% 
% % subplot(2,1,1)
% semilogx(lambda,nn,'k'), hold on
% % semilogx([lambda(2) lambda(end)],[nat nat],'b--'), hold on
% % semilogx(lambda(i),sumPath(i),'r*'), hold off
% % ylim([0 m*n])
% xlabel('\lambda')
% ylabel('accuracy')
% 
% legend({'SINDy path'})


% % figure
% subplot(2,1,2)
% sz=30;
% map = [1 0 0
%     0 0 1];
% % for i = 1:length(lambda)
%     plot(cos(phi),sin(phi),'b--'), hold on
%     plot(cos(phi)*dP,sin(phi)*dP,'r--'), hold on
%     % plot(cos(phi).*XiPath(:,i),sin(phi).*XiPath(:,i),'b','LineWidth',2), hold off
%     scatter(cos(phi).*XiPath(:,i),sin(phi).*XiPath(:,i),sz,XiPath(:,i),'filled'), hold off
%     axis([-1.1 1.1 -1.1 1.1])
%     axis equal
%     axis off
% 
%     colormap(map)
% 
%     txt = ['\lambda = ' num2str(lambda(i),3)];
%     text(0.8,0.8,txt,'FontSize',14);
%     
%     legend({'active terms', 'non-active terms'}, 'Location','Northwest')
%     
%     drawnow
% end