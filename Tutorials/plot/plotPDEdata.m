function plotPDEdata(traj,Ua,name)

dx = traj.dx; 
[Lx,Lt] = size(traj,'uu');

[wT,wX] = meshgrid(traj.tt,dx*(1:Lx));
figure('units','pixels','position',[50 50 750 400])
subplot(1,2,1)
s=surf(wT,wX,traj.uu);
s.EdgeColor = 'none';
xlabel('t'), ylabel('x'), zlabel('u')
title('clean data','interpreter','latex')
subplot(1,2,2)
s=surf(wT,wX,Ua);
s.EdgeColor = 'none';
title('noisy data','interpreter','latex')
xlabel('t'), ylabel('x'), zlabel('u')
sgtitle(sprintf(name),'interpreter','latex')

end