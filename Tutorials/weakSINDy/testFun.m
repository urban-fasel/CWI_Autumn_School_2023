function phi = testFun(k,p,x,t,Hx,Ht)

% dphi(x,p,d) = d/dx^d ((x^2 - 1)^p)
wx = dphi(x,p(1),k(1));   
wt = dphi(t,p(2),k(2));
[wT,wX] = meshgrid(wt,wx);
phi = wX.*wT;

% Define variable conversions: x_underbar to x and t_underbar to t
% chain rule: dphi/dx = dphi/dx_ub * dx_ub/dx
%    dx_ub/dx = 1/dx((x-xk)/Hx) = 1/Hx
% higher derivatives: d^2phi/dx^2 = d^2phi/dx_ub^2 * (dx_ub/dx)^2 + dphi/dx_ub*d^2x_ub/dx^2 
%    d^2x_ub/dx^2 = 0   ->   d^2phi/dx^2 = d^2phi/dx_ub^2 * (dx_ub/dx)^2
S_x = 1/Hx; 
S_t = 1/Ht;  
phi = phi*S_x^k(1)*S_t^k(2);

end

function dphi = dphi(x,p,d)
% dphi = d/dx^d ((x^2 - 1)^p)
dphi = 0;
for k = 0:p
    if 2*(p-k)-d >= 0
        dphi = dphi + nchoosek(p,k)*(-1)^k*factorial(2*(p-k))/factorial(2*(p-k)-d)*x.^(2*(p-k)-d);
    end 
end
end