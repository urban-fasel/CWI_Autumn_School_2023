function dx = lorenz(t,x,param)

sigma = param(1); 
rho = param(2); 
beta = param(3);

dx = [
    sigma*(x(2)-x(1));
    x(1)*(rho-x(3))-x(2);
    x(1)*x(2)-beta*x(3);
    ];