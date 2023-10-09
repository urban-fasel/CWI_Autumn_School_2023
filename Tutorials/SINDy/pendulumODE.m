function xdot = pendulumODE(t,x,u,p) 

% x = [theta; thetaDot]
% u = control

xdot = [-x(2);
        -p.g/p.l*sin(x(1))-p.k/p.m*x(2)+1/(p.m*p.l^2)*u]; % gravity, friction, control