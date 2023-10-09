function plotTF(x,t,phi00,phi01,phi10)

[wT,wX] = meshgrid(t,x);
figure('units','pixels','position',[50 50 900 400])
subplot(1,3,1)
s=surf(wT,wX,phi00);
s.EdgeColor = 'none';
xlabel('t'), ylabel('x'), zlabel('u')
title(sprintf('$(x^2 - 1)^4 (t^2 - 1)^3$'),'interpreter','latex')
subplot(1,3,2)
s=surf(wT,wX,phi01);
s.EdgeColor = 'none';
xlabel('t'), ylabel('x'), zlabel('u')
title(sprintf('$d/dt ( (x^2 - 1)^4 (t^2 - 1)^3 )$'),'interpreter','latex')
subplot(1,3,3)
s=surf(wT,wX,phi10);
s.EdgeColor = 'none';
xlabel('t'), ylabel('x'), zlabel('u')
title(sprintf('$d/dx ( (x^2 - 1)^4 (t^2 - 1)^3 )$'),'interpreter','latex')
sgtitle(sprintf('Test function $\\phi$ and derivatives'),'interpreter','latex')